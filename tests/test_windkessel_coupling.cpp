#include <algorithm>
#include <cmath>
#include <filesystem>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/io/network_parser.hpp"
#include "hemo1d/couplings/register.hpp"
#include "hemo1d/couplings/windkessel_coupling.hpp"
#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/tube_law.hpp"
#include "hemo1d/simulation.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using namespace hemo1d::couplings;
using Catch::Approx;

namespace {

constexpr Real kRho = 1.05;

VesselParameters makeParams() {
    VesselParameters p;
    p.length = 4.0;
    p.A0 = 0.126;
    p.beta = 6.06e5;
    p.alpha = 4.0 / 3.0;
    p.frictionKr = 0.0;
    p.nElements = 64;
    p.polynomialOrder = 1;
    return p;
}

TerminalInterface restInterface(const VesselParameters& p, VesselEnd end = VesselEnd::Distal) {
    TerminalInterface iface;
    iface.end = end;
    iface.params = p;
    iface.trace = {p.A0, 0.0};
    iface.grad = {0.0, 0.0};
    return iface;
}

WindkesselParameters oneCompartment(Real r1, Real r, Real c, Real pOut = 0.0, Real pInit = 0.0) {
    WindkesselParameters wk;
    wk.R1 = r1;
    wk.compartments = {LumpedCompartment{r, c}};
    wk.pOut = pOut;
    wk.pInit = pInit;
    wk.subSteps = 1;
    return wk;
}

// Peak cross-sectional area anywhere in a single-vessel network driven by the
// half-sine inlet flow pulse, with the given windkessel outlet.
Real peakAreaWithOutlet(const WindkesselParameters& wk) {
    const std::filesystem::path csv =
        std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "data" / "inlet_flow_sine.csv";

    VesselParameters vp = makeParams();
    std::vector<Vessel> vessels{Vessel(1, "v", vp)};

    BoundaryConditionSpec inlet;
    inlet.type = BoundaryConditionType::Prescribed;
    inlet.quantity = PrescribedQuantity::FlowRate;
    inlet.csvFile = csv.string();

    BoundaryConditionSpec outlet;
    outlet.type = BoundaryConditionType::NonReflecting; // replaced below

    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, inlet),
        Node(2, "out", {{1, VesselEnd::Distal}}, outlet),
    };
    Network net(FluidProperties{kRho}, std::move(vessels), std::move(nodes));

    SimulationSettings s;
    s.defaultPolynomialOrder = 1;
    s.flux = FluxKind::Hll;
    s.useSlopeLimiter = true;
    Simulation sim(std::move(net), s);
    sim.setCoupling(2, makeWindkesselCoupling(wk));

    const Real dt = 4e-6;
    Real peak = 0.0;
    for (int i = 0; i < 6000; ++i) {
        sim.step(dt);
        if (i % 25 == 0) {
            for (Real a : sim.state().A) {
                REQUIRE(std::isfinite(a));
                REQUIRE(a > 0.0);
                peak = std::max(peak, a);
            }
        }
    }
    for (Real a : sim.state().A) peak = std::max(peak, a);
    return peak;
}

} // namespace

TEST_CASE("WindkesselCoupling leaves an exact rest state unchanged", "[windkessel]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    WindkesselCoupling wk(oneCompartment(/*r1=*/-1.0, /*r=*/5.0e4, /*c=*/1.0e-6));
    const TerminalInterface iface = restInterface(p);

    const SectionState g = wk.solve(iface, model, 0.0, 1e-4);
    CHECK(g.A == Approx(p.A0).margin(1e-9));
    CHECK(g.Q == Approx(0.0).margin(1e-9));

    wk.commit(g, iface, model, 0.0, 1e-4);
    CHECK(wk.compartmentPressures().front() == Approx(0.0).margin(1e-12));
    CHECK(wk.effectiveR1() > 0.0); // matched impedance resolved during solve()
}

TEST_CASE("WindkesselCoupling matched R1 equals rho*c0/A0", "[windkessel]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    WindkesselCoupling wk(oneCompartment(-1.0, 5.0e4, 1.0e-6));
    wk.solve(restInterface(p), model, 0.0, 1e-4);

    CHECK(wk.effectiveR1() == Approx(kRho * model.waveSpeed(p.A0, p) / p.A0));
}

TEST_CASE("WindkesselCoupling commit charges a single RC compartment to Q*R", "[windkessel]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    const Real R = 4.0e4;
    const Real C = 2.0e-6; // R C = 0.08 s
    const Real Q = 3.0;
    const TerminalInterface iface = restInterface(p); // commit only reads iface.end
    const Real dt = 1e-4;

    WindkesselCoupling wk(oneCompartment(/*r1=*/1.0e3, R, C));
    for (int i = 0; i < 4000; ++i) { // t = 0.4 s ~= 5 time constants
        wk.commit(SectionState{p.A0, Q}, iface, model, i * dt, dt);
    }
    CHECK(wk.compartmentPressures().front() == Approx(Q * R).epsilon(0.02));

    WindkesselCoupling half(oneCompartment(1.0e3, R, C));
    const int halfSteps = static_cast<int>(std::lround(R * C * std::log(2.0) / dt));
    for (int i = 0; i < halfSteps; ++i) {
        half.commit(SectionState{p.A0, Q}, iface, model, i * dt, dt);
    }
    CHECK(half.compartmentPressures().front() == Approx(0.5 * Q * R).epsilon(0.03));
}

TEST_CASE("WindkesselCoupling relaxes to p_out once inflow stops", "[windkessel]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    const Real R = 3.0e4;
    const Real C = 2.0e-6;
    const Real pOutVal = 150.0;
    WindkesselCoupling wk(oneCompartment(1.0e3, R, C, pOutVal, /*pInit=*/900.0));

    const TerminalInterface iface = restInterface(p);
    const Real dt = 1e-4;
    for (int i = 0; i < 5000; ++i) { // no inflow: Q = 0
        wk.commit(SectionState{p.A0, 0.0}, iface, model, i * dt, dt);
    }
    CHECK(wk.compartmentPressures().front() == Approx(pOutVal).epsilon(0.02));
}

TEST_CASE("WindkesselCoupling solve enforces the RCR interface relation", "[windkessel]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    const Real R1 = 6.0e3;
    const Real P1 = 400.0;
    WindkesselCoupling wk(oneCompartment(R1, 5.0e4, 1.0e-6, /*pOut=*/0.0, /*pInit=*/P1));

    const Real aStar = p.A0 * 1.03;
    const Real pAstar = model.pressure(aStar, p);
    const Real dp = law.pressureDerivative(aStar, p);

    SECTION("arbitrary trace: the linearized windkessel row holds exactly") {
        TerminalInterface iface = restInterface(p);
        iface.trace = {aStar, 0.02};
        const SectionState g = wk.solve(iface, model, 0.0, 1e-4);
        // dp*A - R1*Q = P1 - p(A*) + dp*A*  <=>  p_lin(A) - R1*Q = P1
        CHECK(dp * g.A - R1 * g.Q == Approx(P1 - pAstar + dp * aStar).margin(1e-6));
    }

    SECTION("an RCR fixed point is returned unchanged") {
        const Real qStar = (pAstar - P1) / R1; // p(A*) - P1 = R1 Q*  -> steady RCR
        TerminalInterface iface = restInterface(p);
        iface.trace = {aStar, qStar}; // grad = 0
        const SectionState g = wk.solve(iface, model, 0.0, 1e-4);
        CHECK(g.A == Approx(aStar).margin(1e-9));
        CHECK(g.Q == Approx(qStar).margin(1e-9));
    }
}

TEST_CASE("WindkesselCoupling: a resistive bed reflects a forward pulse (positive)", "[windkessel]") {
    const VesselParameters p = makeParams();
    const Real matchedR1 = kRho * std::sqrt(p.beta / (2.0 * kRho * std::sqrt(p.A0))) / p.A0;

    const Real matchedPeak = peakAreaWithOutlet(oneCompartment(-1.0, 2.0e3, 5.0e-6));
    const Real occludedPeak =
        peakAreaWithOutlet(oneCompartment(20.0 * matchedR1, 5.0e5, 1.0e-7));

    CHECK(occludedPeak > matchedPeak);  // the resistive bed piles area up
    CHECK(occludedPeak > p.A0 * 1.02);  // a clear positive (compression) reflection
    CHECK(matchedPeak < p.A0 * 1.05);   // matched bed: little reflection
}

TEST_CASE("BoundaryResolver populates TerminalInterface::rho for couplings", "[windkessel]") {
    // A custom coupling (here a CallbackCoupling) must be able to form its own
    // wave speeds -- so the resolver hands it the fluid density.
    registerBuiltinCouplings(); // the example JSON now uses "external"/"windkessel"
    Network net =
        io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "single_vessel_windkessel.json");
    SimulationSettings s;
    s.flux = FluxKind::Hll;
    Simulation sim(std::move(net), s);

    Real seenRho = -1.0;
    sim.setCouplingCallback(2, [&](const TerminalInterface& iface, Real, Real) -> SectionState {
        seenRho = iface.rho;
        return iface.trace; // rigid-ish passthrough, keeps the run finite
    });
    sim.step(4e-6);
    CHECK(seenRho == Approx(sim.network().fluid().density));
    CHECK(seenRho > 0.0);
}

TEST_CASE("Simulation drives the single-vessel windkessel example and stays stable", "[windkessel]") {
    registerBuiltinCouplings();
    Network net =
        io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "single_vessel_windkessel.json");

    SimulationSettings s;
    s.defaultPolynomialOrder = 1;
    s.flux = FluxKind::Hll;
    s.useSlopeLimiter = true;
    Simulation sim(std::move(net), s);
    sim.addProbe("mid", 1, 2.0);

    const Real dt = 4e-6;
    sim.run(4000 * dt, dt, /*recordEvery=*/10);

    for (Real a : sim.state().A) {
        REQUIRE(std::isfinite(a));
        REQUIRE(a > 0.0);
    }
    for (Real q : sim.state().Q) REQUIRE(std::isfinite(q));
    CHECK(sim.probeSamples("mid").size() == 400);
}