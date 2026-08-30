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

def write_network(name, h, p, log=False):
    network_path = os.path.join(c.GENERATED_DIR, name)
    c.write_network_json(network_path, h, p, c.PULSE_CSV_NAME)
    if log:
        print(f"Wrote {network_path}")
    return network_path

def load_network(network_path):
    network = hemo1d.load_network(network_path)
    print(f"Loaded: {network.vessel_count} vessels, {network.node_count} nodes")
    return network

def run_case(study, case_name, network_path, p, dt, h, lengths):
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
        study=study, case_name=case_name, polynomial_order=p,
        dt=dt, h=h, cfl_margin=margin,
        target_time=c.TARGET_TIME,
        record_every=record_every, steps=steps,
        num_snapshots=snapshots, lengths=lengths
    )
    print(f"  [{study}/{case_name}] h={h} p={p} dt={dt:.4e} steps={steps} "
          f"cfl={margin:.2%} snapshots={snapshots} -> {elapsed:.1f}s")


def run_bifurcation_cases():
    for p in c.P_LIST:
        for h in c.H_LIST:
            name = c.format_case_name(p, h)
            path = os.path.join(c.GENERATED_DIR, name + ".json")
            c.write_bifurcation_network_json(path, h, p, c.PULSE_CSV_NAME)
            c.check_cfl(c.DT, h, p, f"{c.BIFURCATION_SUBDIR}/{name}")
            run_case(c.BIFURCATION_SUBDIR, name, path, p, c.DT, h, 
                     [c.BIFURCATION_PARENT_LENGTH, c.BIFURCATION_DAUGHTER_LENGTH, c.BIFURCATION_DAUGHTER_LENGTH])

def run_single_cases():
    for p in c.P_LIST:    
        for h in c.H_LIST:
            name = c.format_case_name(p, h)
            network_path = write_network(name+".json", h, p)
            run_case(c.SPATIAL_SUBDIR, name, network_path, p, c.DT, h, [c.LENGTH,])
            

def main():
    write_pulse()
    run_single_cases()  
    run_bifurcation_cases()
    

# def main():
#     h, p = 0.125, 1
#     name = f"Y-p{p}_h{h:g}".replace(".", "p")     # p1_h0p125
#     network_path = os.path.join(c.GENERATED_DIR, name + ".json")
#     c.write_bifurcation_network_json(network_path, h, p, c.PULSE_CSV_NAME)
#     run_case("bifurcation", name, network_path, p, c.DT, h)


if __name__ == "__main__":
    main()