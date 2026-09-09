import os
import argparse
import glob
import math
import numpy as np

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
    parser.add_argument("--dir", default = "output_windkessel")
    parser.add_argument("--quantity", default="flow_rate", choices=sorted(QUANTITY_LABELS))
    parser.add_argument("--cols", type=int, default=3)
    parser.add_argument("--t-min", type=float, default=None, help="Only plot t >= this (e.g. to skip the initial transient).")
    args = parser.parse_args()

    probes_dir = os.path.join(THIS_DIR, args.dir, "probes")
    csv_paths = sorted(glob.glob(os.path.join(probes_dir, "v*.csv")))
    if not csv_paths:
        raise AttributeError(f"No probe CSV found in ", probes_dir)
    min_global = float("inf")
    max_global = -1*float("inf")
    for csv_path in csv_paths:
        df = pd.read_csv(csv_path)
        min_global = min(df[args.quantity].min(), min_global)
        max_global = max(df[args.quantity].max(), max_global)
    print(f"Range {args.dir} {args.quantity}:[{min_global},{max_global}]")
    
if __name__ == "__main__":
    main()