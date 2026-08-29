import os
import math
import json

RHO = 1.055
VISCOSITY = 0.045
A0 = 0.126
BETA = 606060.0
ALPHA = 4.0 / 3.0
LENGTH = 4.0

KINEMATIC_VISCOSITY = VISCOSITY / RHO
FRICTION_KR = 8.0 * math.pi * KINEMATIC_VISCOSITY

C0 = math.sqrt(BETA * math.sqrt(A0) / (2.0 * RHO * A0))

PULSE_AMPLITUDE = 5.0
PULSE_HALF_WIDTH = 4.0e-4
PULSE_CENTER = 1.0e-3
PULSE_WINDOW = 2.0 * PULSE_CENTER
PULSE_CSV_DT = 2.0e-7

TARGET_TIME = 3.0e-3
NUM_SNAPSHOTS = 100
BASE_DT_INTERVAL = TARGET_TIME / NUM_SNAPSHOTS

CONVERGENCE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR      = os.path.join(CONVERGENCE_DIR, "output")
GENERATED_DIR   = os.path.join(OUTPUT_DIR, "_generated")
PULSE_CSV_NAME  = "inlet_pulse.csv"

def pulse_value(t):
    if abs(t-PULSE_CENTER) > PULSE_HALF_WIDTH:
        return 0.0
    return PULSE_AMPLITUDE * math.sin( \
        math.pi * (t - PULSE_CENTER + PULSE_HALF_WIDTH) / (2 * PULSE_HALF_WIDTH))**2

def write_pulse_csv(path):
    n = round(PULSE_WINDOW / PULSE_CSV_DT)
    with open(path, "w") as f:
        f.write("time, flow_rate\n")
        for i in range(n+1):
            t = i*PULSE_CSV_DT
            f.write(f"{t:.9f}, {pulse_value(t):.8f}\n")

def element_size(n_elements):
    return LENGTH / n_elements

def cfl_max_dt(h, p):
    return h / ((2*p+1) * C0)

def check_cfl(dt, h, p, label):
    # dt margin
    safe_dt = cfl_max_dt(h, p)
    assert dt < safe_dt, f"{label}: dt: {dt:.4e} is bigger than cfl: {safe_dt:.4e}"
    return dt / safe_dt

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

def build_network(n_elements, polynomial_order, pulse_csv_name):
     # A single vessel: prescribed flow-rate inlet, non-reflecting outlet
    return {
        "_description": (
            f"Convergence-study single vessel "
            f"(n_elements={n_elements}, polynomial_order={polynomial_order})"
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
                "n_elements": n_elements,
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


def write_network_json(path, n_elements, polynomial_order, pulse_csv_name):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(build_network(n_elements, polynomial_order, pulse_csv_name), f, indent=2)