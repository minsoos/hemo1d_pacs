import os
import sys

import hemo1d
import matplotlib.pyplot as plt


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(THIS_DIR)

QUANTITIES = [
    ("area", "Area [cm^2]"),
    ("flow_rate", "Flow rate [cm^3/s]"),
    ("pressure", "Pressure (P - P0) [dyn/cm^2]"),
    ("velocity", "Velocity [cm/s]")
]


##### EXAMPLE WITH RAMP AREA #####
# EXAMPLE_PARAMS = {
#     "network": "simple_bifurcation.json",
#     "target_time": 3e-3,
#     "output": "out_ramp_area"
# }

##### EXAMPLE WITH SINE FLOW RATE #####
EXAMPLE_PARAMS = {
    "network": "simple_bifurcation_sine.json",
    "target_time": 0.3,
    "output": "out_sine_flow"
}



def vessel_of(probe_name):
    return probe_name.rsplit("_", 1)[0] if "_" in probe_name else probe_name


def plot_probes(sim, probe_names, output_dir, title_prefix):
    vessels = {}
    for name in probe_names:
        vessels.setdefault(vessel_of(name), []).append(name)

    for vessel, names in vessels.items():
        fig, axes = plt.subplots(len(QUANTITIES), 1, sharex=True, figsize=(8, 8))

        for ax, (attr, ylabel) in zip(axes, QUANTITIES):
            for name in names:
                samples = sim.probe_samples(name)

                time = [s.time for s in samples]
                values = [getattr(s, attr) for s in samples]

                ax.plot(time, values, label=name)
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.3)

        axes[0].legend(loc="best", fontsize="small")
        axes[-1].set_xlabel("time [s]")
        fig.suptitle(f"{title_prefix} -- vessel {vessel}")
        fig.tight_layout()

        plot_path = os.path.join(output_dir, f"probes_{vessel}.png")
        fig.savefig(plot_path, dpi=150)
        print(f"Plot saved to {plot_path}")

    plt.show()


def main():
    output_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(THIS_DIR, EXAMPLE_PARAMS["output"])
    os.makedirs(output_dir, exist_ok=True)

    ##### LOAD NETWORK #####
    network_path = os.path.join(THIS_DIR, "networks", EXAMPLE_PARAMS["network"])
    network = hemo1d.load_network(network_path)
    print(f"Loaded network: {network.vessel_count} vessels, {network.node_count} nodes.")

    ##### SET SIM SETTINGS #####
    settings = hemo1d.SimulationSettings()
    settings.default_polynomial_order = 1
    settings.flux = hemo1d.FluxKind.HLL
    settings.use_slope_limiter = True

    ##### CREATE SIMULATION #####
    sim = hemo1d.Simulation(network, settings)

    ##### ADD PROBES #####
    sim.add_probe("omega1_inlet", vessel_id=1, z=0.05)
    sim.add_probe("omega1_mid", vessel_id=1, z=0.5)
    sim.add_probe("omega2_mid", vessel_id=2, z=0.5)
    sim.add_probe("omega3_mid", vessel_id=3, z=1.0)
    sim.add_probe("omega2_outlet", vessel_id=2, z=0.95)
    sim.add_probe("omega3_outlet", vessel_id=3, z=1.95)

    ##### SET SIM PARAMS #####
    dt = 5e-6
    target_time = EXAMPLE_PARAMS["target_time"]
    vtk_dir = os.path.join(output_dir, "vtk")

    ##### RUN SIMULATION #####
    print(f"Running to t={target_time:g}s with dt={dt:g}s...")
    sim.run(target_time, dt, record_every=2, vtk_directory=vtk_dir, vtk_every=20)
    print(f"Done. Final simulation time: {sim.time:g}s")

    ##### WRITE OUTPUTS #####
    sim.write_probes_csv(output_dir)
    print(f"Probe CSVs written to {output_dir}.")
    print(f"VTK time series written to {vtk_dir} (open {os.path.join(vtk_dir, 'hemo1d.pvd')} in ParaView)")

    ##### PLOT PROBES #####
    probe_names = ["omega1_inlet", "omega1_mid", "omega2_mid", "omega2_outlet", "omega3_mid", "omega3_outlet"]
    plot_probes(sim, probe_names, output_dir, title_prefix="Hemo1D probes: simple bifurcation")


if __name__ == "__main__":
    main()