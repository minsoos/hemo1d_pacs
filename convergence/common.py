import os
import math
import json
import numpy as np


# %% Physic parameters ------------------------------------------
RHO = 1.055
VISCOSITY = 0.045
A0 = 0.126
BETA = 606060.0
ALPHA = 4.0 / 3.0

KINEMATIC_VISCOSITY = VISCOSITY / RHO
FRICTION_KR = 8.0 * math.pi * KINEMATIC_VISCOSITY

C0 = math.sqrt(BETA * math.sqrt(A0) / (2.0 * RHO * A0))

# %% Pulses -----------------------------------------------------

# sin^2 pulse
PULSE_AMPLITUDE = 5.0
PULSE_HALF_WIDTH = 4.0e-4
PULSE_CENTER = 1.0e-3
PULSE_WINDOW = 2.0 * PULSE_CENTER
PULSE_CSV_DT = 2.0e-7

def pulse_value(t):
    if abs(t-PULSE_CENTER) > PULSE_HALF_WIDTH:
        return 0.0
    return PULSE_AMPLITUDE * math.sin( \
        math.pi * (t - PULSE_CENTER + PULSE_HALF_WIDTH) / (2 * PULSE_HALF_WIDTH))**2


# gaussian pulse
PULSE_AMPLITUDE = 5.0
PULSE_SIGMA     = 2.0e-4
PULSE_CENTER    = 5.0 * PULSE_SIGMA
PULSE_WINDOW    = 2.0 * PULSE_CENTER
PULSE_CSV_DT    = 2.0e-7

def pulse_value(t):
    return PULSE_AMPLITUDE * math.exp(-0.5 * ((t - PULSE_CENTER) / PULSE_SIGMA) ** 2)

def write_pulse_csv(path):
    n = round(PULSE_WINDOW / PULSE_CSV_DT)
    with open(path, "w") as f:
        f.write("time, flow_rate\n")
        for i in range(n+1):
            t = i*PULSE_CSV_DT
            f.write(f"{t:.9f}, {pulse_value(t):.8f}\n")

# %% Time and snapshots
TARGET_TIME = 3.0e-3
NUM_SNAPSHOTS = 100
BASE_DT_INTERVAL = TARGET_TIME / NUM_SNAPSHOTS

def total_steps_for(dt):
    n = TARGET_TIME / dt
    n_int = round(n)
    assert abs(n-n_int) < 1e-6 * max(n, 1.0), f"dt={dt:.4e} does not evenly divide TARGET_TIME"
    return n_int

def record_steps_for(dt):
    k = BASE_DT_INTERVAL / dt
    k_int = round(k)
    assert abs(k-k_int) < 1e-6 * max(k, 1.0), f"dt={dt:.4e} does not evenly divide BASE_DT_INTERVAL"
    return k_int

# %% Directories
CONVERGENCE_DIR = os.path.dirname(os.path.abspath(__file__))
SPATIAL_SUBDIR     = "spatial"
BIFURCATION_SUBDIR = "bifurcation"
OUTPUT_DIR      = os.path.join(CONVERGENCE_DIR, "output")
GENERATED_DIR   = os.path.join(OUTPUT_DIR, "_generated")
PULSE_CSV_NAME  = "inlet_pulse.csv"

def format_case_name(p, h):
    return f"p{p}_h{h:g}".replace(".", "p")





# %% Convergence refinement parameters
P_LIST   = [1, 2]
LENGTH = 4.0
H_LIST = tuple(1/np.pow(2,np.arange(3,9))) ##
DT       = BASE_DT_INTERVAL / 384


FIELD_CSV_COLUMNS = ["snapshot_index", "snapshot_time", "element_index",
                     "vessel_id", "z", "area", "flow_rate"]




def element_size(n_elements):
    return LENGTH / n_elements

def cfl_max_dt(h, p):
    return h / ((2*p+1) * C0)

def check_cfl(dt, h, p, label):
    # dt margin
    safe_dt = cfl_max_dt(h, p)
    assert dt < safe_dt, f"{label}: dt: {dt:.4e} is bigger than cfl: {safe_dt:.4e}"
    return dt / safe_dt


# %% Single network ----------------------------------
def build_network(h, polynomial_order, pulse_csv_name):
     # A single vessel: prescribed flow-rate inlet, non-reflecting outlet
    n = round(LENGTH/h)
    return {
        "_description": (
            f"Convergence-study single vessel "
            f"(n_elements={n}, polynomial_order={polynomial_order})"
        ),
        "fluid": {"density": RHO, "viscosity": VISCOSITY},
        "vessels": [
            {
                "id": 1,
                "name": "vessel",
                "length": LENGTH,
                "A0": A0,
                "beta": BETA,
                "alpha": ALPHA,
                "n_elements": n,
                "polynomial_order": polynomial_order,
            }
        ],
        "nodes": [
            {
                "id": 1, "name": "inlet",
                "connections": [{"vessel": 1, "end": "proximal"}],
                "boundary_condition": {
                    "type": "prescribed",
                    "quantity": "flow_rate",
                    "csv_file": pulse_csv_name,
                },
            },
            {
                "id": 2, "name": "outlet",
                "connections": [{"vessel": 1, "end": "distal"}],
                "boundary_condition": {"type": "non_reflecting"},
            },
        ],
    }

def write_network_json(path, h, polynomial_order, pulse_csv_name):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(build_network(h, polynomial_order, pulse_csv_name), f, indent=2)


# %% Y-shape network
BIFURCATION_PARENT_LENGTH   = 1.0
BIFURCATION_DAUGHTER_LENGTH = 2.0

def bifurcation_elements(h):
    return round(BIFURCATION_PARENT_LENGTH / h), round(BIFURCATION_DAUGHTER_LENGTH / h)

def build_bifurcation_network(h, polynomial_order, pulse_csv_name):
    n_par, n_dau = bifurcation_elements(h)

    def vessel(vid, name, length, n):
        return {"id": vid, "name": name, "length": length, "A0": A0,
                "beta": BETA, "alpha": ALPHA, "n_elements": n,
                "polynomial_order": polynomial_order}

    return {
        "_description": f"Convergence-study symmetric bifurcation (h={h}, p={polynomial_order})",
        "fluid": {"density": RHO, "viscosity": VISCOSITY},
        "vessels": [
            vessel(1, "parent",     BIFURCATION_PARENT_LENGTH,   n_par),
            vessel(2, "daughter_a", BIFURCATION_DAUGHTER_LENGTH, n_dau),
            vessel(3, "daughter_b", BIFURCATION_DAUGHTER_LENGTH, n_dau),
        ],
        "nodes": [
            {"id": 1, "name": "inlet",
             "connections": [{"vessel": 1, "end": "proximal"}],
             "boundary_condition": {"type": "prescribed", "quantity": "flow_rate",
                                    "csv_file": pulse_csv_name}},
            {"id": 2, "name": "junction",
             "connections": [{"vessel": 1, "end": "distal"},
                             {"vessel": 2, "end": "proximal"},
                             {"vessel": 3, "end": "proximal"}],
            "bifurcation_angles_rad": [0.7853981634, 0.7853981634]},
            {"id": 3, "name": "outlet_a",
             "connections": [{"vessel": 2, "end": "distal"}],
             "boundary_condition": {"type": "non_reflecting"}},
            {"id": 4, "name": "outlet_b",
             "connections": [{"vessel": 3, "end": "distal"}],
             "boundary_condition": {"type": "non_reflecting"}},
        ],
    }


def write_bifurcation_network_json(path, h, polynomial_order, pulse_csv_name):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(build_bifurcation_network(h, polynomial_order, pulse_csv_name), f, indent=2)


# %% Manifest

def case_dir(study, case_name):
    d = os.path.join(OUTPUT_DIR, study, case_name)
    os.makedirs(d, exist_ok=True)
    return d

def write_manifest(path, **fields):
    with open(path, "w") as f:
        json.dump(fields, f, indent=2)

def read_manifest(path):
    with open(path, "r") as f:
        return json.load(f)

def main():
    import hemo1d
    network_path = os.path.join(GENERATED_DIR, "Y-p1_n0008.json")
    write_bifurcation_network_json(network_path, 0.125, 1, PULSE_CSV_NAME)
    network = hemo1d.load_network(network_path)
    print(f"Loaded: {network.vessel_count} vessels, {network.node_count} nodes")
    return network



if __name__ == "__main__":
    main()