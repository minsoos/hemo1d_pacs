import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)

import build_network as bn  # noqa: E402  -- VESSELS / INFLOWS / OUTLETS / periods

MMHG = 1333.22  # dyn/cm^2 per mmHg

COMPLETE_DIR = os.path.join(THIS_DIR, "output_windkessel", "probes")
INCOMPLETE_DIR = os.path.join(THIS_DIR, "output_windkessel_incomplete", "probes")

INLETS = set(bn.INFLOWS)          # {1, 5, 12}
OUTLETS = set(bn.OUTLETS)         # {3, 8, 10, 15, 17, 18}
COMMUNICATING = {4, 11, 16}       # L-PComm, R-PComm, ACoA -- the collateral ring

GROUP_COLOR = {
    "inlet": "#4C72B0",
    "outlet": "#C44E52",
    "communicating": "#55A868",
    "other": "#9C9C9C",
}


def group_of(vid):
    if vid in INLETS:
        return "inlet"
    if vid in OUTLETS:
        return "outlet"
    if vid in COMMUNICATING:
        return "communicating"
    return "other"


def short_label(vid):
    """Compact anatomical tag: L-MCA, R-PCA II, ACoA, BA, ..."""
    n = bn.VESSELS[vid][0].lower()
    side = "L" if n.startswith("left") else "R" if n.startswith("right") else ""
    core = ("PComm" if "posterior communicating" in n else
            "ACoA" if "anterior communicating" in n else
            "PCA" if "posterior cerebral" in n else
            "MCA" if "middle cerebral" in n else
            "ACA" if "anterior cerebral" in n else
            "ICA" if "internal carotid" in n else
            "BA" if "basilar" in n else bn.VESSELS[vid][0])
    seg = " II" if n.endswith(" ii") else " I" if n.endswith(" i") else ""
    return f"{side}-{core}{seg}".lstrip("-")


def load_run(probes_dir):
    """{vessel_id -> DataFrame(time, area, flow_rate, pressure, velocity)}."""
    if not os.path.isdir(probes_dir):
        raise SystemExit(
            f"{probes_dir}\n  not found -- run the simulations first "
            f"(see this script's header)."
        )
    runs = {}
    for fn in sorted(os.listdir(probes_dir)):
        if not (fn.startswith("v") and fn.endswith(".csv")):
            continue
        vid = int(fn[1:].split("_", 1)[0])
        runs[vid] = pd.read_csv(os.path.join(probes_dir, fn))
    if not runs:
        raise SystemExit(f"no v*.csv probe files in {probes_dir}")
    return runs


def window_stats(df, t0, t1):
    d = df[(df["time"] >= t0) & (df["time"] <= t1)]
    q = d["flow_rate"].to_numpy()
    p = d["pressure"].to_numpy() / MMHG
    return dict(
        q_mean=float(q.mean()),
        q_absmean=float(np.abs(q).mean()),
        q_min=float(q.min()),
        q_max=float(q.max()),
        q_amp=float(q.max() - q.min()),
        p_mean=float(p.mean()),
        p_amp=float(p.max() - p.min()),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t-min", type=float, default=bn.CARDIAC_PERIOD,
                    help="start of the averaging window [s] (default: one cardiac period, "
                         "i.e. discard the first cycle).")
    ap.add_argument("--t-max", type=float, default=None,
                    help="end of the averaging window [s] (default: last common sample).")
    ap.add_argument("--vessels", default=None,
                    help="comma-separated vessel ids for the waveform figure "
                         "(default: the 6 outlets + the communicating arteries still present).")
    ap.add_argument("--outdir", default=THIS_DIR)
    args = ap.parse_args()

    complete = load_run(COMPLETE_DIR)
    incomplete = load_run(INCOMPLETE_DIR)

    removed = sorted(set(complete) - set(incomplete))
    common = sorted(set(complete) & set(incomplete))

    # Common averaging window.
    last = min(float(complete[v]["time"].iloc[-1]) for v in common)
    last = min(last, min(float(incomplete[v]["time"].iloc[-1]) for v in common))
    t0 = args.t_min
    t1 = args.t_max if args.t_max is not None else last
    n_cycles = (t1 - t0) / bn.CARDIAC_PERIOD
    print(f"removed vessel(s): {[f'{v} ({bn.VESSELS[v][0]})' for v in removed] or 'none'}")
    print(f"averaging window: [{t0:.2f}, {t1:.2f}] s  (~{n_cycles:.2f} cardiac cycles)\n")

    # ---- per-vessel table -------------------------------------------------
    rows = []
    for vid in common:
        c = window_stats(complete[vid], t0, t1)
        i = window_stats(incomplete[vid], t0, t1)
        dq = i["q_absmean"] - c["q_absmean"]
        pct = 100.0 * dq / c["q_absmean"] if c["q_absmean"] > 1e-9 else np.nan
        rows.append(dict(
            vid=vid, vessel=bn.VESSELS[vid][0], tag=short_label(vid), group=group_of(vid),
            q_mean_complete=c["q_mean"], q_mean_incomplete=i["q_mean"],
            q_absmean_complete=c["q_absmean"], q_absmean_incomplete=i["q_absmean"],
            d_q_absmean=dq, pct_q_absmean=pct,
            q_amp_complete=c["q_amp"], q_amp_incomplete=i["q_amp"],
            p_mean_mmHg_complete=c["p_mean"], p_mean_mmHg_incomplete=i["p_mean"],
            d_p_mean_mmHg=i["p_mean"] - c["p_mean"],
            p_amp_mmHg_complete=c["p_amp"], p_amp_mmHg_incomplete=i["p_amp"],
        ))
    table = pd.DataFrame(rows).sort_values("d_q_absmean")
    csv_path = os.path.join(args.outdir, "incomplete_effect_summary.csv")
    table.to_csv(csv_path, index=False, float_format="%.5g")

    def show(sub, title):
        if sub.empty:
            return
        print(title)
        for _, r in sub.iterrows():
            print(f"  {r.tag:<9} q_mean {r.q_mean_complete:8.3f} -> {r.q_mean_incomplete:8.3f} "
                  f"cm3/s   |q| {r.pct_q_absmean:+6.1f}%   "
                  f"p_mean {r.p_mean_mmHg_complete:6.1f} -> {r.p_mean_mmHg_incomplete:6.1f} mmHg")
        print()

    show(table[table.group == "inlet"], "Inlets")
    show(table[table.group == "communicating"], "Communicating arteries (collateral ring)")
    show(table[table.group == "outlet"], "Territorial outlets (perfusion)")
    movers = table[(table.group == "other")].reindex(
        table[(table.group == "other")].d_q_absmean.abs().sort_values(ascending=False).index).head(6)
    show(movers, "Largest changes among the remaining vessels")

    # ---- conservation / total perfusion --------------------------------
    def total_abs(run):
        return sum(window_stats(run[v], t0, t1)["q_absmean"] for v in OUTLETS if v in run)

    perf_c, perf_i = total_abs(complete), total_abs(incomplete)
    print(f"total outlet flow (sum of 6 territories):  "
          f"{perf_c:.2f} -> {perf_i:.2f} cm3/s  ({100 * (perf_i - perf_c) / perf_c:+.1f}%)")
    print("  (inlets are prescribed velocities, so total supply is ~fixed; the change "
          "is redistribution, plus any compliance/steady-state drift.)\n")

    # ---- Fig A: per-vessel change in mean |flow| -----------------------
    fa = table.sort_values("d_q_absmean")
    fig, ax = plt.subplots(figsize=(7.5, 0.32 * len(fa) + 1.2))
    ax.barh(fa.tag, fa.d_q_absmean, color=[GROUP_COLOR[g] for g in fa.group])
    ax.axvline(0.0, color="k", lw=0.8)
    ax.set_xlabel(r"$\Delta$ mean $|$flow$|$  (incomplete $-$ complete)  [cm$^3$/s]")
    ax.set_title(f"Per-vessel flow change -- {', '.join(short_label(v) for v in removed)} removed",
                 fontsize=11)
    ax.grid(axis="x", alpha=0.3)
    handles = [plt.Rectangle((0, 0), 1, 1, color=GROUP_COLOR[g]) for g in GROUP_COLOR]
    ax.legend(handles, GROUP_COLOR.keys(), fontsize=8, loc="lower right")
    fig.tight_layout()
    pa = os.path.join(args.outdir, "incomplete_effect_flow.png")
    fig.savefig(pa, dpi=150)

    # ---- Fig B: the six territorial outlets ---------------------------
    outs = [v for v in bn.OUTLETS if v in common]
    by_vid = table.set_index("vid")
    xc = np.arange(len(outs))
    qc = [by_vid.loc[v, "q_absmean_complete"] for v in outs]
    qi = [by_vid.loc[v, "q_absmean_incomplete"] for v in outs]
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    ax.bar(xc - 0.2, qc, 0.4, label="complete", color="#4C72B0")
    ax.bar(xc + 0.2, qi, 0.4, label="incomplete", color="#C44E52")
    for k, v in enumerate(outs):
        pct = 100.0 * (qi[k] - qc[k]) / qc[k]
        ax.annotate(f"{pct:+.1f}%", (xc[k], max(qc[k], qi[k])), textcoords="offset points",
                    xytext=(0, 3), ha="center", fontsize=8)
    ax.set_xticks(xc)
    ax.set_xticklabels([short_label(v) for v in outs])
    ax.set_ylabel(r"mean $|$flow$|$ at vessel midpoint  [cm$^3$/s]")
    ax.set_title("Territorial perfusion: complete vs incomplete Circle of Willis", fontsize=11)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.15)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    pb = os.path.join(args.outdir, "incomplete_effect_outlets.png")
    fig.savefig(pb, dpi=150)

    # ---- Fig C: waveform overlays -----------------------------------
    if args.vessels:
        sel = [int(x) for x in args.vessels.split(",")]
    else:
        sel = [v for v in bn.OUTLETS] + sorted(COMMUNICATING & set(common))
    sel = [v for v in sel if v in common]
    ncol = 3
    nrow = int(np.ceil(len(sel) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.5 * nrow), squeeze=False)
    for ax, vid in zip(axes.flat, sel):
        dc = complete[vid]
        di = incomplete[vid]
        ax.plot(dc["time"], dc["flow_rate"], lw=0.9, color="#4C72B0", label="complete")
        ax.plot(di["time"], di["flow_rate"], lw=0.9, color="#C44E52", ls="--", label="incomplete")
        ax.axvspan(t0, t1, color="k", alpha=0.05)
        ax.set_title(f"{short_label(vid)}  ({vid})", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)
    for ax in axes.flat[len(sel):]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=8)
    fig.supxlabel("time [s]", fontsize=9)
    fig.supylabel(r"flow rate [cm$^3$/s]", fontsize=9)
    fig.suptitle("Flow waveforms at vessel midpoints -- shaded = averaging window", fontsize=11)
    fig.tight_layout()
    pc = os.path.join(args.outdir, "incomplete_effect_waveforms.png")
    fig.savefig(pc, dpi=150)

    print("wrote:")
    for p in (csv_path, pa, pb, pc):
        print("  " + os.path.relpath(p, args.outdir))


if __name__ == "__main__":
    main()
