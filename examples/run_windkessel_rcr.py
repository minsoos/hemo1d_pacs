import math
import os
import sys

import hemo1d


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
NETWORK = os.path.join(THIS_DIR, "networks", "windkessel_rcr.json")

OUTLET_NODE = 2
R2, C = 5.0e4, 4.0e-6          # the RCR params in windkessel_rcr.json
DT, T_END, PULSE_END = 2.0e-5, 1.0, 0.25
MMHG = 1333.22                 # dyn/cm^2 per mmHg


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(THIS_DIR, "out_windkessel_rcr")
    os.makedirs(out_dir, exist_ok=True)

    settings = hemo1d.SimulationSettings()
    settings.default_polynomial_order = 1
    settings.flux = hemo1d.FluxKind.HLL
    settings.use_slope_limiter = True
    sim = hemo1d.Simulation(hemo1d.load_network(NETWORK), settings)
    sim.add_probe("inlet", vessel_id=1, z=0.1)
    sim.add_probe("outlet", vessel_id=1, z=9.7)

    # Step in chunks so the 0D compartment pressure P_c can be sampled over time
    # (probe histories are filled by run(); P_c is read with coupling_state()).
    n_steps, rec = int(round(T_END / DT)), 100
    pc = []
    for k in range(0, n_steps, rec):
        sim.run(min((k + rec) * DT, T_END), DT, record_every=rec)
        pc.append((sim.time, sim.coupling_state(OUTLET_NODE)[0]))

    inl = sim.probe_samples("inlet")
    o = sim.probe_samples("outlet")

    for name, samples in (("inlet", inl), ("outlet", o)):
        with open(os.path.join(out_dir, f"{name}.csv"), "w") as fh:
            fh.write("time,area,flow_rate,pressure\n")
            for s in samples:
                fh.write(f"{s.time},{s.area},{s.flow_rate},{s.pressure}\n")
    with open(os.path.join(out_dir, "compartment.csv"), "w") as fh:
        fh.write("time,p_compartment\n")
        for t, p in pc:
            fh.write(f"{t},{p}\n")
    print(f"\nCSVs -> {out_dir}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib not available -- skipping the plot)")
        return

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(8, 8))
    ax[0].plot([s.time for s in inl], [s.flow_rate for s in inl], "--", color="#8A928C",
               label="inlet flow (prescribed pulse)")
    ax[0].plot([s.time for s in o], [s.flow_rate for s in o], "-", color="#3F7A4E",
               label="outlet flow (vessel -> R-C bed)")
    ax[0].set_ylabel("flow  Q  [cm^3/s]")
    ax[1].plot([s.time for s in o], [s.pressure for s in o], color="#2E6E8C")
    ax[1].set_ylabel("outlet pressure\nP - P0  [dyn/cm^2]")
    ax[2].plot([t for t, _ in pc], [p for _, p in pc], color="#8A6D1E")
    ax[2].set_ylabel("0D compartment\npressure P_c  [dyn/cm^2]")
    ax[2].set_xlabel("time  [s]")
    for a in ax:
        a.axvspan(0, PULSE_END, color="k", alpha=0.06, label="systole (inflow on)")
        a.axvline(PULSE_END, color="k", lw=0.8, ls="--")
        a.grid(alpha=0.3)
    ax[0].legend(loc="upper right", fontsize=9)
    fig.suptitle("Three-element Windkessel outlet")
    fig.tight_layout()
    path = os.path.join(out_dir, "windkessel_rcr.png")
    fig.savefig(path, dpi=140)
    print(f"plot -> {path}")
    plt.show()


if __name__ == "__main__":
    main()
