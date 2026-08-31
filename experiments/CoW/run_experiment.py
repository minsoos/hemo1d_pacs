import os
import sys
import argparse
import hemo1d
import time
import math

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)

import build_network as bn

GENERATED_DIR = os.path.join(THIS_DIR, "generated")
OUTPUT_DIR = os.path.join(THIS_DIR, "output")

def slug(name):
    return name.replace(" ", "_").replace("-", "_")

def check_state(sim):
    snap = sim.field_snapshot()
    bad = [s for s in snap
           if not (math.isfinite(s.area) and s.area > 0.0 and math.isfinite(s.flow_rate))]
    if bad:
        s = bad[0]
        raise RuntimeError(
            f"invalid state at t={sim.time:.6e}: {len(bad)}/{len(snap)} DOFs bad. "
            f"First at vessel {s.vessel_id}, z={s.z:.4f}, A={s.area:.4f}, Q={s.flow_rate:.4f}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-time", type=float, default=bn.SIMULATION_DURATION,
                         help="Total simulated time in seconds (default: full 2-cycle run).")
    parser.add_argument("--cfl-number", type=float, default=0.8,
                         help="Safety factor for the dynamic CFL-based dt.")
    parser.add_argument("--record-every", type=int, default=300,
                         help="Record every probe every N steps (~2ms cadence at the measured dt).")
    parser.add_argument("--vtk-every", type=int, default=3000,
                         help="Write a VTK snapshot every N steps (~20ms cadence). 0 disables VTK.")
    parser.add_argument("--polynomial-order", type=int, default=bn.POLYNOMIAL_ORDER)
    args = parser.parse_args()

    network_path = os.path.join(GENERATED_DIR, "cow_network.json")
    if not os.path.exists(network_path):
        raise SystemExit(f"{network_path} not found. Run build_network.py first")
    network = hemo1d.load_network(network_path)

    settings = hemo1d.SimulationSettings()
    settings.default_polynomial_order = args.polynomial_order
    settings.flux = hemo1d.FluxKind.HLL
    settings.use_slope_limiter = True

    sim = hemo1d.Simulation(network, settings)

    for vid, (name, length, a0, beta) in sorted(bn.VESSELS.items()):
        sim.add_probe(f"v{vid:02d}_{slug(name)}", vid, length/2.0)

    probes_dir = os.path.join(OUTPUT_DIR, "probes")
    os.makedirs(probes_dir, exist_ok=True)
    vtk_dir = os.path.join(OUTPUT_DIR, "vtk") if args.vtk_every>0 else ""

    if vtk_dir:
        os.makedirs(vtk_dir, exist_ok=True)

    initial_dt = sim.cfl_time_step(args.cfl_number)
    print(f"Simulation begins: A priori CFL:", initial_dt)
    t0 = time.time()
    sim.run(args.target_time, -1.0, args.record_every, vtk_dir, args.vtk_every, args.cfl_number)
    elapsed = time.time() - t0
    print(f"Simulation finished up to: {sim.time:.6f}s. Task duration: {elapsed:.1f}s")

    check_state(sim)

    sim.write_probes_csv(probes_dir)
    print(f"Wrote {len(bn.VESSELS)} probe CSVs to {probes_dir}")
    if vtk_dir:
        print(f"Wrote VTK series to {vtk_dir}")


if __name__ == "__main__":
    main()