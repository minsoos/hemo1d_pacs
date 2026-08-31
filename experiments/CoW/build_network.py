import os
import math
import csv
import json

# %% Directories
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
INFLOWS_DIR = os.path.join(THIS_DIR, "inflows")
GENERATED_DIR = os.path.join(THIS_DIR, "generated")

# %% Physical parameters
DENSITY = 1.055
VISCOSITY = 0.045
ALPHA = 4.0 / 3.0
FRICTION_KR = 8.0 * math.pi * (VISCOSITY / DENSITY)

H = 0.125
DT = 1e-5
T = 1
POLYNOMIAL_ORDER = 2

# Table 1: vessels -> (name, length, A0, beta)
VESSELS = {
    1: ("Basilar artery", 1.2800, 0.0740, 0.1451e7),
    2: ("Left posterior cerebral artery I", 1.3212, 0.0347, 0.0497e7),
    3: ("Left posterior cerebral artery II", 3.0909, 0.0338, 0.0490e7),
    4: ("Left posterior communicating artery", 1.1713, 0.0256, 0.0426e7),
    5: ("Left internal carotid artery I", 4.2887, 0.1252, 0.1887e7),
    6: ("Left internal carotid artery II", 0.4203, 0.0756, 0.0733e7),
    7: ("Left anterior cerebral artery I", 3.0665, 0.0288, 0.0452e7),
    8: ("Left middle cerebral artery", 4.4854, 0.0580, 0.0642e7),
    9: ("Right posterior cerebral artery I", 1.0012, 0.0295, 0.0458e7),
    10: ("Right posterior cerebral artery II", 2.5496, 0.0331, 0.0485e7),
    11: ("Right posterior communicating artery", 1.8936, 0.0301, 0.0463e7),
    12: ("Right internal carotid artery I", 4.8650, 0.1260, 0.1893e7),
    13: ("Right internal carotid artery II", 0.6961, 0.0800, 0.0754e7),
    14: ("Right anterior cerebral artery I", 2.1127, 0.0435, 0.0556e7),
    15: ("Right middle cerebral artery", 3.3325, 0.0623, 0.0666e7),
    16: ("Anterior communicating artery", 0.4632, 0.0214, 0.0390e7),
    17: ("Left anterior cerebral artery II", 4.4019, 0.0281, 0.0447e7),
    18: ("Right anterior cerebral artery II", 3.9169, 0.0323, 0.0479e7),
}


# Table 2: bifurcation -> ((id1, end1), (id2, end2), (id3, end3), (theta2, theta3)).
BIFURCATIONS = {
    "I": ((1, "distal"), (2, "proximal"), (9, "proximal"), (0.7748, 0.7771)),
    "II": ((2, "distal"), (3, "proximal"), (4, "proximal"), (0.7054, 0.9872)),
    "III": ((5, "distal"), (4, "distal"), (6, "proximal"), (1.2182, 0.5778)),
    "IV": ((6, "distal"), (7, "proximal"), (8, "proximal"), (1.9317, 0.1917)),
    "V": ((7, "distal"), (16, "proximal"), (17, "proximal"), (2.4837, 0.2554)),
    "VI": ((9, "distal"), (10, "proximal"), (11, "proximal"), (0.8728, 1.1379)),
    "VII": ((12, "distal"), (13, "proximal"), (11, "distal"), (0.6948, 1.3204)),
    "VIII": ((13, "distal"), (14, "proximal"), (15, "proximal"), (1.0024, 0.9820)),
    "IX": ((14, "distal"), (16, "distal"), (18, "proximal"), (0.8910, 0.9184)),
}

# Inlet inflows
INFLOWS = {
    1: "BAS.csv",
    5: "L-ICA_I.csv",
    12: "R-ICA_I.csv",
}

# Outlet vessels
OUTLETS = [3, 8, 10, 15, 17, 18]

SIMULATION_DURATION = 2.0  # s, matches the thesis' own run length.
CARDIAC_PERIOD = 1.0  # s, the inflow CSVs cover exactly one period each.


def make_periodic_csv(in_path, out_path, period, n_cycles):
    times = []
    values = []
    with open(in_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            times.append(float(row[0]))
            values.append(float(row[1]))
    t0 = times[0]
    shifted = [t - t0 for t in times]
    assert abs(shifted[-1] - period) < 1e-6, f"{in_path}: last sample is not one period after the first"
    assert abs(values[0] - values[-1]) < 1e-9, f"{in_path}: data is not periodic (first != last sample)"

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time", "value"])
        for cycle in range(n_cycles):
            offset = cycle * period
            start_idx = 1 if cycle > 0 else 0 # Skip the first value
            # So it doesn't repeat with the last of previous period
            for i in range(start_idx, len(shifted)):
                writer.writerow([f"{shifted[i] + offset:.10e}", f"{values[i]:.10e}"])


def make_constant_csv(path, value, end_time):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time", "value"])
        writer.writerow([0.0, value])
        writer.writerow([end_time, value])


def build_network_json():
    n_cycles = math.ceil(SIMULATION_DURATION / CARDIAC_PERIOD)

    # Build vessels
    vessels_json = []
    for vid, (name, length, a0, beta) in sorted(VESSELS.items()):
        n_elements = math.ceil(length / H)
        vessels_json.append(
            {
                "id": vid,
                "name": name.replace(" ", "_"),
                "length": length,
                "A0": a0,
                "beta": beta,
                "alpha": ALPHA,
                "friction_kr": FRICTION_KR,
                "n_elements": n_elements,
            }
        )

    # Build nodes
    nodes_json = []
    next_node_id = 1

    for vid, csv_name in sorted(INFLOWS.items()):
        periodic_name = csv_name.strip(".csv") + "_periodic.csv"
        nodes_json.append(
            {
                "id": next_node_id,
                "name": f"inlet_{VESSELS[vid][0].replace(' ', '_')}",
                "connections": [{"vessel": vid, "end": "proximal"}],
                "boundary_condition": {
                    "type": "prescribed",
                    "quantity": "velocity",
                    "csv_file": f"periodic_inflows/{periodic_name}",
                },
            }
        )
        next_node_id += 1

    for vid in OUTLETS:
        nodes_json.append(
            {
                "id": next_node_id,
                "name": f"outlet_{VESSELS[vid][0].replace(' ', '_')}",
                "connections": [{"vessel": vid, "end": "distal"}],
                "boundary_condition": {
                    "type": "prescribed",
                    "quantity": "pressure",
                    "csv_file": "outlet_pressure_zero.csv",
                },
            }
        )
        next_node_id += 1

    for bif_id, (c1, c2, c3, angles) in BIFURCATIONS.items():
        nodes_json.append(
            {
                "id": next_node_id,
                "name": f"bifurcation_{bif_id}",
                "connections": [
                    {"vessel": c1[0], "end": c1[1]},
                    {"vessel": c2[0], "end": c2[1]},
                    {"vessel": c3[0], "end": c3[1]},
                ],
                "bifurcation_angles_rad": list(angles),
            }
        )
        next_node_id += 1

    network = {
        "_description": (
            "Complete Circle of Willis network (18 vessels, 9 trifurcations, 3 "
            "prescribed-velocity inlets, 6 prescribed-pressure outlets), geometry "
        ),
        "fluid": {"density": DENSITY, "viscosity": VISCOSITY},
        "vessels": vessels_json,
        "nodes": nodes_json,
    }
    return network, n_cycles


def main():
    os.makedirs(GENERATED_DIR, exist_ok=True)
    periodic_dir = os.path.join(GENERATED_DIR, "periodic_inflows")
    os.makedirs(periodic_dir, exist_ok=True)

    network, n_cycles = build_network_json()

    for vid, csv_name in INFLOWS.items():
        in_path = os.path.join(INFLOWS_DIR, csv_name)
        periodic_name = csv_name.strip(".csv") + "_periodic.csv"
        out_path = os.path.join(periodic_dir, periodic_name)
        make_periodic_csv(in_path, out_path, CARDIAC_PERIOD, n_cycles)
        print(f"Wrote {out_path} ({n_cycles} cycles, {n_cycles * CARDIAC_PERIOD:.1f}s)")

    outlet_csv = os.path.join(GENERATED_DIR, "outlet_pressure_zero.csv")
    make_constant_csv(outlet_csv, 0.0, SIMULATION_DURATION + 1.0)
    print(f"Wrote {outlet_csv}")

    network_path = os.path.join(GENERATED_DIR, "cow_network.json")
    with open(network_path, "w") as f:
        json.dump(network, f, indent=2)
    print(f"Wrote {network_path}")

    total_elements = sum(v["n_elements"] for v in network["vessels"])
    print(f"\n{len(network['vessels'])} vessels, {len(network['nodes'])} nodes, "
          f"{total_elements} elements total (target h={H} cm, p={POLYNOMIAL_ORDER}).")

    
if __name__ == "__main__":
    main()