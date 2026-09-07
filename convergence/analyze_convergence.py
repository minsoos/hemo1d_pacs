import csv
import math
import os
import sys
import argparse
 
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as c


def load_case(case_dir, requested_time=None):
    manifest = c.read_manifest(os.path.join(case_dir, "manifest.json"))

    num_snapshots = int(manifest["num_snapshots"])
    final_time = float(manifest["target_time"])    
    

    if requested_time is None:
        selected_index = num_snapshots
    else:
        if requested_time <= 0.0 or requested_time > final_time:
            raise ValueError(f"load_case: Requested time must be in (0,{final_time}]")
        snapshot_interval  = final_time / num_snapshots
        selected_index = int(
            math.floor(requested_time / snapshot_interval + 0.5)
        )
        selected_index = max(1, min(num_snapshots, selected_index))


    by_vessel = {}
    with open(os.path.join(case_dir, "field_snapshots.csv"), newline="") as f:
        for row in csv.DictReader(f):
            if int(row["snapshot_index"]) != selected_index:
                continue

            
            vid = int(row["vessel_id"])
            e = int(row["element_index"])
            by_vessel.setdefault(vid, {}).setdefault(e, []).append(
                (float(row["z"]), float(row["area"]), float(row["flow_rate"]))
            )
    if not by_vessel:
        raise RuntimeError(f"Snapshot {selected_index} was not found")
    return manifest, by_vessel

class PiecewiseField:

    def __init__(self, n_elements, length, elem_dict, quantity_index):
        self.n = n_elements
        self.h = length / n_elements
        self.length = length
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

def l2_error(fields_a, fields_b, h, degree):
    xi, w = np.polynomial.legendre.leggauss(degree + 1)
    per_vessel = {}
    for vid in sorted(fields_a):
        fa, fb = fields_a[vid], fields_b[vid]
        assert abs(fa.length - fb.length) < 1e-12, "l2_error: length mismatch"
        n_quad = round(fa.length / h)
        total = 0.0
        for k in range(n_quad):
            mid, half = k * h + 0.5 * h, 0.5 * h
            for xi_i, w_i in zip(xi, w):
                z = mid + half * xi_i
                d = fa(z) - fb(z)
                total += w_i * d * d * half
        per_vessel[vid] = math.sqrt(total)
    return per_vessel

def field_for_vessel(elems):
    n = len(elems)
    length = max(z for pts in elems.values() for (z, _, _) in pts)
    return PiecewiseField(n, length, elems, _QUANTITY_INDEX)

def _field(h, p, requested_time=None):
    _, by_vessel = load_case(os.path.join(c.OUTPUT_DIR, STUDY, c.format_case_name(p,h)), requested_time)
    return {vid: field_for_vessel(elems)
                      for vid, elems in by_vessel.items()}

class Richardson:
    def __init__(self, fine_field, coarse_field, order):
        assert fine_field.length==coarse_field.length, "Richardson: Fields must be same length"
        self.fine = fine_field
        self.length = fine_field.length
        self.coarse = coarse_field
        self.factor = 1.0 / (2.0 ** order - 1.0)

    def __call__(self, z):
        uf = self.fine(z)
        uc = self.coarse(z)
        return uf + (uf-uc) * self.factor
    

def reference_field(h, p, order, requested_time=None):
    fields_fine = _field(h/2, p, requested_time)
    fields_coarse = _field(h, p, requested_time)
    return {id: Richardson(fields_fine[id], fields_coarse[id], order) for id in
        fields_fine}

def loglog_order(hs, errs):
    return float(np.polyfit(np.log(np.array(hs[-3:])), np.log(np.array(errs[-3:])), 1)[0])

def global_norm(per_vessel):
    return math.sqrt(sum(e * e for e in per_vessel.values()))

def analyze_convergence(p, hs=(1/8,1/16,1/32), requested_time=None, expected_order=None, richardson=True):
    if expected_order is None:
        expected_order = p+1
    if richardson:
        ref = reference_field(hs[-1]/2, p, order=expected_order, requested_time=requested_time)
    else:
        ref = _field(hs[-1]/2.0, p, requested_time=requested_time)
    per_h = []

    for h in hs:
        field = _field(h, p, requested_time=requested_time)
        per_h.append(l2_error(field, ref, hs[-1]/(2)**2, p))
    errs  = [global_norm(d) for d in per_h]
    return list(hs), errs, per_h, loglog_order(hs, errs)


def normalize_list(x):
    x = np.array(x)
    return x/x[0]


def main():
    global STUDY, QUANTITY, _QUANTITY_INDEX

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--study",
        choices=(
            c.SPATIAL_SUBDIR,
            c.BIFURCATION_SUBDIR,
            c.WINDKESSEL_SUBDIR,
        ),
        required=True,
    )
    parser.add_argument(
        "--quantity",
        choices=("area", "flow_rate"),
        default="flow_rate",
    )
    parser.add_argument(
        "--time",
        type=float,
        default=None,
        help="Physical snapshot time; default is the final snapshot.",
    )
    parser.add_argument(
        "--expected_order",
        type=int,
        default=None,
        help="Expected order, to pass to Richardson.",
    )
    parser.add_argument(
            "--not_richardson_use",
            action="store_true",
            help="Expected order, to pass to Richardson.",
        )

    parser.add_argument(
            "--extended_hs",
            action="store_true",
            help="Use extended refinemenet.",
        )
    args = parser.parse_args()

    STUDY = args.study
    QUANTITY = args.quantity
    _QUANTITY_INDEX = {
        "area": 0,
        "flow_rate": 1,
    }[QUANTITY]
    REQ_TIME = args.time
    EXPECTED_ORDER = args.expected_order
    RICHARDSON_USE = not args.not_richardson_use

    h_list = c.H_LIST
    if args.extended_hs:
        h_list = c.H_LIST_EXTENDED

    if RICHARDSON_USE:
        hs = h_list[:-2]
    else:
        hs = h_list[:-1]
    print("hs:", hs)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.set_xticks(hs)
    ax.set_xticklabels([f"{h:g}" for h in hs])
    ax.xaxis.set_minor_formatter(NullFormatter())

    print(f"{'p':>2} {'h':>10} {'error':>14} {'order':>8}")
    results = {}
    for p in c.P_LIST:
        expected_order = EXPECTED_ORDER
        if EXPECTED_ORDER is None:
            expected_order = p+1
        hs_out, errs, per_h, order = analyze_convergence(p, hs, requested_time=REQ_TIME,
                                                         expected_order=expected_order,
                                                         richardson=RICHARDSON_USE)
        results[p] = (hs_out, errs, per_h, order)

        for i, (h, e) in enumerate(zip(hs, errs)):
            ratio = errs[i-1] / e if i else float("nan")
            order_obs = np.log2(ratio)
            print(f"{p:>2} {h:>10.5f} {e:>14.6e} {order_obs:>8.2f}")
        for vid, e in sorted(per_h[-1].items()):
            print(f"      vaso {vid}: {e:.6e}")
        print(f"   p={p}: observed order {order:.3f} (expected {expected_order})\n")
        ax.loglog(hs, errs, "o-", label=f"p={p} (order {order:.2f})")

    for p, q, style in [(1, 2, "--"), (2, 3, ":")]:
        hs_out, errs, _, _ = results[p]
        h_ref, e_ref = hs_out[-1], errs[-1]
        ax.loglog(hs_out, e_ref * (np.array(hs_out) / h_ref) ** q,
                  style, color="0.5", lw=1, label=f"$h^{q}$")
        print([f"{e:.3e}" for e in errs])
        print([c.format_case_name(p, h) for h in hs])

    ax.set_xlabel("h [cm]")
    ax.set_ylabel(f"L2 error in {QUANTITY}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()

    time_tag = "final" if REQ_TIME is None else f"{REQ_TIME:g}"
    out = os.path.join(c.OUTPUT_DIR, f"{STUDY}_order_{QUANTITY}_T{time_tag}.png")
    fig.savefig(out, dpi=150)
    print("wrote", out)

if __name__ == "__main__":
    main()

