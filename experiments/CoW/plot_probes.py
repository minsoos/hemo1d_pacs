import os
import argparse
import glob
import math

import matplotlib.pyplot as plt
import pandas as pd

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

QUANTITY_LABELS = {
    "area": "Area [cm^2]",
    "flow_rate": "Flow rate [cm^3/s]",
    "pressure": "Pressure [g/(cm*s^2)]",
    "velocity": "Velocity [cm/s]",
}

def vessel_label(csv_path):
    stem = os.path.splitext(os.path.basename(csv_path))[0]
    vid, name = stem.split("_", 1)
    return vid.lstrip("v").lstrip("0") or "0", name.replace("_", " ")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probes-dir", default=os.path.join(THIS_DIR, "output", "probes"))
    parser.add_argument("--quantity", default="flow_rate", choices=sorted(QUANTITY_LABELS))
    parser.add_argument("--output", default=os.path.join(THIS_DIR, "output", "probes_grid.png"))
    parser.add_argument("--cols", type=int, default=3)
    parser.add_argument("--t-min", type=float, default=None, help="Only plot t >= this (e.g. to skip the initial transient).")
    args = parser.parse_args()

    csv_paths = sorted(glob.glob(os.path.join(args.probes_dir, "v*.csv")))
    if not csv_paths:
        raise AttributeError(f"No probe CSV found in ", args.probes_dir)

    n = len(csv_paths)
    cols = args.cols
    rows = math.ceil(n/cols)

    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 2.6 * rows), sharex=True)
    axes = axes.flatten()

    for ax, csv_path in zip(axes, csv_paths):
        vid, name = vessel_label(csv_path)
        df = pd.read_csv(csv_path)
        if args.t_min is not None:
            df = df[df["time"] >= args.t_min]
        ax.plot(df["time"], df[args.quantity], linewidth=0.8)
        ax.set_title(f"{vid}: {name}", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

    for ax in axes[n:]:
        ax.axis("off")

    for ax in axes[max(0, n - cols):n]:
        ax.set_xlabel("time [s]", fontsize=8)

    fig.suptitle(f"CoW probes -- {QUANTITY_LABELS[args.quantity]} at each vessel midpoint", fontsize=12)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")

    
if __name__ == "__main__":
    main()