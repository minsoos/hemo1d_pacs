import os
import sys
import hemo1d
import time
import csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import common as c

def write_pulse(log=False):
    os.makedirs(c.GENERATED_DIR, exist_ok=True)
    pulse_csv_path = os.path.join(c.GENERATED_DIR, c.PULSE_CSV_NAME)
    c.write_pulse_csv(pulse_csv_path)
    if log:
        print(f"Wrote inlet pulse CSV: {pulse_csv_path}")

def write_network(name, n, p, log=False):
    network_path = os.path.join(c.GENERATED_DIR, name)
    c.write_network_json(network_path, n, p, c.PULSE_CSV_NAME)
    if log:
        print(f"Wrote {network_path}")
    return network_path

def load_network(network_path):
    network = hemo1d.load_network(network_path)
    print(f"Loaded: {network.vessel_count} vessels, {network.node_count} nodes")
    return network

def run_case(study, case_name, network_path, n, p, dt, h):
    out_dir = c.case_dir(study, case_name)
    field_csv_path = os.path.join(out_dir, "field_snapshots.csv")
    manifest_path = os.path.join(out_dir, "manifest.json")

    margin = c.check_cfl(dt, h, p, out_dir)
    steps = c.total_steps_for(dt)
    record_every = c.record_steps_for(dt)

    network = hemo1d.load_network(network_path)

    settings = hemo1d.SimulationSettings()
    settings.default_polynomial_order = p
    settings.flux = hemo1d.FluxKind.HLL
    settings.use_slope_limiter = False

    sim = hemo1d.Simulation(network, settings)

    t0 = time.time()
    snapshots = 0
    with open(field_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(c.FIELD_CSV_COLUMNS)
        for step_index in range(1, steps+1):
            sim.step(dt)
            if step_index % record_every == 0:
                snapshots += 1
                snapshot_time = sim.time
                for s in sim.field_snapshot():
                     writer.writerow([snapshots, snapshot_time, s.element_index,
                                     s.vessel_id, s.z, s.area, s.flow_rate])
    elapsed = time.time() - t0

    c.write_manifest(
        manifest_path,
        study=study, case_name=case_name,
        n_elements=n, polynomial_order=p,
        dt=dt, h=h, cfl_margin=margin,
        target_time=c.TARGET_TIME,
        record_every=record_every, steps=steps,
        num_snapshots=snapshots
    )
    print(f"  [{study}/{case_name}] n={n} p={p} dt={dt:.4e} steps={steps} "
          f"cfl={margin:.2%} snapshots={snapshots} -> {elapsed:.1f}s")




def main():
    write_pulse()
    dt = c.SPATIAL_DT
    
    for p in c.SPATIAL_P_LIST:    
        for n in c.SPATIAL_N_LIST:
            h = c.element_size(n)
            network_name = f"p{p}_n{n:04d}.json"
            network_path = write_network(network_name, n, p)
            run_case("spatial", f"p{p}_n{n:04d}", network_path, n, p, dt, h)
        
    

if __name__ == "__main__":
    main()