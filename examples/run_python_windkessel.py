import math
import os
import sys

import hemo1d

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
NETWORK = os.path.join(THIS_DIR, "networks", "single_vessel_windkessel.json")
OUTLET_NODE = 2
DT, TARGET_T = 4e-6, 0.30
PROBES = [("inlet", 0.05), ("mid", 2.0), ("outlet", 3.95)]

# The single_vessel_windkessel.json outlet, as plain numbers so both the
# built-in and the Python coupling get exactly the same parameters.
WK = dict(r1=-1.0,                                   # < 0 -> matched impedance rho*c0/A0
          compartments=[(5.0e4, 8.0e-7), (3.0e4, 5.0e-7)],
          p_out=0.0, p_init=0.0, sub_steps=2)


# --------------------------------------------------------------------------
# The coupling, entirely in Python. Mirrors hemo1d::couplings::WindkesselCoupling for the
# linear-elastic tube law  P - P0 = beta * (sqrt(A) - sqrt(A0)) / A0.
# --------------------------------------------------------------------------
class PythonWindkessel:
    def __init__(self, compartments, r1=-1.0, p_out=0.0, p_init=0.0, sub_steps=1):
        self.RC = [tuple(c) for c in compartments]   # [(R_k, C_k), ...], interface-first
        self.r1_cfg = float(r1)
        self.p_out = float(p_out)
        self.sub = int(sub_steps)
        self.P = [float(p_init)] * len(self.RC)       # compartment pressures = the 0D state
        self.r1 = None

    # linear-elastic tube law helpers
    @staticmethod
    def _p(A, A0, beta):     return beta * (math.sqrt(A) - math.sqrt(A0)) / A0     # P - P0
    @staticmethod
    def _dpdA(A, A0, beta):  return beta / (2.0 * A0 * math.sqrt(A))               # dP/dA
    @staticmethod
    def _c(A, A0, beta, rho):                                                      # wave speed
        return math.sqrt(beta * math.sqrt(A) / (2.0 * rho * A0))

    def _ensure_r1(self, iface):
        if self.r1 is None:
            self.r1 = (self.r1_cfg if self.r1_cfg > 0.0
                       else iface.rho * self._c(iface.A0, iface.A0, iface.beta, iface.rho) / iface.A0)

    def solve(self, iface, t, dt):
        self._ensure_r1(iface)
        A0, beta, alpha, rho = iface.A0, iface.beta, iface.alpha, iface.rho
        As, Qs = iface.A, iface.Q
        u = Qs / As

        # outgoing (distal) left eigenvector of the flux Jacobian
        c = self._c(As, A0, beta, rho)
        c_alpha = math.sqrt(c * c + alpha * (alpha - 1.0) * u * u)
        lout = (c_alpha - alpha * u, 1.0)                        # (l0, l1)

        # forward-Euler compatibility prediction  U* - dt (H dU/dz + B)
        hA = iface.dQdz
        hQ = (c * c - alpha * u * u) * iface.dAdz + 2.0 * alpha * u * iface.dQdz
        bQ = iface.friction_kr * u
        ccA, ccQ = As - dt * hA, Qs - dt * (hQ + bQ)
        rhs1 = lout[0] * ccA + lout[1] * ccQ                     # compatibility row

        # windkessel closure row, linearized about As:
        #   p(A) - R1 Q = P_1   ~=>   dp*A - R1*Q = P_1 - p0 + dp*As
        dp, p0 = self._dpdA(As, A0, beta), self._p(As, A0, beta)
        m00, m01, rhs0 = dp, -self.r1, self.P[0] - p0 + dp * As

        # solve  [m00 m01; l0 l1] (A, Q) = (rhs0, rhs1)
        det = m00 * lout[1] - m01 * lout[0]
        A_out = (rhs0 * lout[1] - m01 * rhs1) / det
        Q_out = (m00 * rhs1 - rhs0 * lout[0]) / det
        return (A_out, Q_out)

    def commit(self, resolved, iface, t, dt):
        q_interface = resolved[1]                                # flow into compartment 1 (distal)
        h, n = dt / self.sub, len(self.RC)
        for _ in range(self.sub):
            qout = [0.0] * n
            for k in range(n):
                nxt = self.P[k + 1] if k + 1 < n else self.p_out
                qout[k] = (self.P[k] - nxt) / self.RC[k][0]
            for k in range(n):
                qin = q_interface if k == 0 else qout[k - 1]
                self.P[k] += h / self.RC[k][1] * (qin - qout[k])


# --------------------------------------------------------------------------
def make_sim():
    settings = hemo1d.SimulationSettings()
    settings.default_polynomial_order = 1
    settings.flux = hemo1d.FluxKind.HLL
    settings.use_slope_limiter = True
    sim = hemo1d.Simulation(hemo1d.load_network(NETWORK), settings)
    for name, z in PROBES:
        sim.add_probe(name, vessel_id=1, z=z)
    return sim


def run(sim):
    sim.run(TARGET_T, DT, record_every=20)
    return {name: sim.probe_samples(name) for name, _ in PROBES}


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(THIS_DIR, "out_python_windkessel")
    os.makedirs(out_dir, exist_ok=True)

    # A: built-in C++ WindkesselCoupling
    sim_a = make_sim()
    sim_a.set_windkessel_outlet(OUTLET_NODE, **WK)
    probes_a = run(sim_a)
    p1_a = sim_a.coupling_state(OUTLET_NODE)          # C++ compartment pressures

    # B: the Python coupling above, same parameters
    sim_b = make_sim()
    wk = PythonWindkessel(**WK)
    sim_b.set_coupling_callback(OUTLET_NODE, wk.solve, wk.commit)
    probes_b = run(sim_b)
    p1_b = list(wk.P)                                 # Python compartment pressures

    # --- agreement check ---------------------------------------------------
    def max_abs_rel(a, b, attr):
        m_abs = m_rel = 0.0
        for sa, sb in zip(a, b):
            d = abs(getattr(sa, attr) - getattr(sb, attr))
            m_abs = max(m_abs, d)
            m_rel = max(m_rel, d / (1.0 + abs(getattr(sa, attr))))
        return m_abs, m_rel

    print("built-in C++ WindkesselCoupling  vs  pure-Python PythonWindkessel")
    print(f"  same params: r1=matched, compartments={WK['compartments']}, p_out=0, sub_steps=2")
    for attr, unit in (("pressure", "dyn/cm^2"), ("flow_rate", "cm^3/s"), ("area", "cm^2")):
        a_abs, a_rel = max_abs_rel(probes_a["outlet"], probes_b["outlet"], attr)
        print(f"  outlet {attr:10s}: max |diff| = {a_abs:.3e} {unit:9s}  (rel {a_rel:.2e})")
    print(f"  compartment pressures  C++ = {[round(x, 4) for x in p1_a]}")
    print(f"                      Python = {[round(x, 4) for x in p1_b]}")

    # --- plot ------------------------------------------------------------
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not available -- skipping the plot)")
        return

    fig, ax = plt.subplots(3, 1, sharex=True, figsize=(9, 9))
    for attr, a, ylabel in (("pressure", ax[0], "outlet P - P0 [dyn/cm^2]"),
                            ("flow_rate", ax[1], "outlet Q [cm^3/s]")):
        sa, sb = probes_a["outlet"], probes_b["outlet"]
        a.plot([s.time for s in sa], [getattr(s, attr) for s in sa], "-", lw=2.2,
               color="#2E6E8C", label="built-in C++ WindkesselCoupling")
        a.plot([s.time for s in sb], [getattr(s, attr) for s in sb], "--", lw=2.2,
               color="#A8583B", label="pure-Python PythonWindkessel")
        a.set_ylabel(ylabel)
        a.grid(True, alpha=0.3)
    ax[0].legend(loc="best", fontsize=9)
    # difference on its own axis
    sa, sb = probes_a["outlet"], probes_b["outlet"]
    ax[2].plot([s.time for s in sa],
               [sa[i].pressure - sb[i].pressure for i in range(len(sa))], color="#8A6D1E")
    ax[2].set_ylabel("P difference  (C++ - Python)\n[dyn/cm^2]")
    ax[2].set_xlabel("time [s]")
    ax[2].grid(True, alpha=0.3)
    fig.suptitle("Same 0D Windkessel, two implementations: C++ built-in vs pure-Python callback")
    fig.tight_layout()
    path = os.path.join(out_dir, "python_windkessel.png")
    fig.savefig(path, dpi=140)
    print(f"\nPlot saved to {path}")
    plt.show()


if __name__ == "__main__":
    main()
