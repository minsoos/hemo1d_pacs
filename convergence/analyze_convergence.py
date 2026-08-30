import csv
import math
import os
import sys
 
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as c

QUANTITY = "flow_rate"
_QUANTITY_INDEX = {"area": 0, "flow_rate": 1}[QUANTITY]

GAUSS_ORDER = 6

def load_case(case_dir):
    manifest = c.read_manifest(os.path.join(case_dir, "manifest.json"))
    last = manifest["num_snapshots"]
    elems = {}
    with open(os.path.join(case_dir, "field_snapshots.csv"), newline="") as f:
        for row in csv.DictReader(f):
            if int(row["snapshot_index"]) != last:
                continue

            e = int(row["element_index"])
            elems.setdefault(e, []).append(
                (float(row["z"]), float(row["area"]), float(row["flow_rate"]))
            )
    return manifest, elems

class PiecewiseField:

    def __init__(self, n_elements, length, elem_dict, quantity_index):
        self.n = n_elements
        self.h = length / n_elements
        sorted_eids = sorted(elem_dict.keys())
        assert len(sorted_eids) == n_elements, (
            f"expected {n_elements} elements, got {len(sorted_eids)}"
        )

        self.polys = [None] * n_elements
        for local_i, eid in enumerate(sorted_eids):
            pts = sorted(elem_dict[eid], key=lambda p: p[0]) # left and right
            z = np.array([p[0] for p in pts])
            v = np.array([p[quantity_index+1] for p in pts])
            degree = len(z)-1
            self.polys[local_i] = np.poly1d(np.polyfit(z, v, degree))

    def __call__(self, z):
        idx = min(self.n-1, max(0, int(z/self.h)))
        return self.polys[idx](z)

def l2_error(field_a, field_b, n_quad_elements, length, degree):
    xi, w = np.polynomial.legendre.leggauss(degree+1)
    h = length / n_quad_elements
    total = 0.0
    for k in range(n_quad_elements):
        a = k*h
        mid, half = a + 0.5*h, 0.5*h
        for xi_i, w_i in zip(xi,w):
            z = mid + half*xi_i
            d = field_a(z) - field_b(z)
            total += w_i * d * d * half
    return math.sqrt(total)

def _field(n,p):
    _, elems = load_case(os.path.join(c.OUTPUT_DIR, "spatial", f"p{p}_n{n:04d}"))
    return PiecewiseField(n, c.LENGTH, elems, _QUANTITY_INDEX)

class Richardson:
    def __init__(self, fine_field, coarse_field, order):
        self.fine = fine_field
        self.coarse = coarse_field
        self.factor = 1.0 / (2.0 ** order - 1.0)

    def __call__(self, z):
        uf = self.fine(z)
        uc = self.coarse(z)
        return uf + (uf-uc) * self.factor

def reference_field(n, p, order):
    field_fine = _field(n*2, p)
    field_coarse = _field(n, p)
    return Richardson(field_fine, field_coarse, order)

def loglog_order(hs, errs):
    return float(np.polyfit(np.log(np.array(hs[-3:])), np.log(np.array(errs[-3:])), 1)[0])

def analyze_spatial(p, levels=(8,16,32)):
    ref = reference_field(levels[-1]*2, p, order=p+1)
    hs, errs = [],[]

    for n in levels:
        field = _field(n, p)
        hs.append(c.LENGTH / n)
        errs.append(l2_error(field, ref, levels[-1]*(2)**2, c.LENGTH, p))
    return hs, errs, loglog_order(hs, errs)

def normalize_list(x):
    x = np.array(x)
    return x/x[0]


def main():
    levels = c.SPATIAL_N_LIST[:-2]
    hs = c.LENGTH/np.array(levels)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.set_xticks(hs)
    ax.set_xticklabels([f"{h:g}" for h in hs])
    ax.xaxis.set_minor_formatter(NullFormatter())
    for q, style in [(2, "--"), (3, ":")]:
        ax.loglog(hs, normalize_list(np.pow(hs,q)), style, label=f"x^{q}")
    print(f"{'p':>2} {'h':>10} {'error':>14} {'ratio':>8}")
    for p in c.SPATIAL_P_LIST:
        hs, errs, order = analyze_spatial(p, levels=levels)
        for i, (h, e) in enumerate(zip(hs, errs)):
            ratio = errs[i-1] / e if i else float("nan")
            print(f"{p:>2} {h:>10.5f} {e:>14.6e} {ratio:>8.2f}")
        print(f"   p={p}: observed order {order:.3f} (expected {p+1})\n")
        ax.loglog(hs, normalize_list(errs), "o-", label=f"p={p} (order {order:.2f})")


    ax.set_xlabel("h [cm]")
    ax.set_ylabel(f"L2 error in {QUANTITY}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = os.path.join(c.OUTPUT_DIR, f"spatial_order_{QUANTITY}.png")
    fig.savefig(out, dpi=150)
    print("wrote", out)

if __name__ == "__main__":
    main()




def build_network_fields_over_time(frames, vessel_specs, quantity_index):
    return {
        idx: {
            vid: PiecewiseField(n, length, by_vessel[vid], quantity_index)
                for vid, (n, length) in vessel_specs.items()
        } for idx, by_vessel in frames.items()
    }

