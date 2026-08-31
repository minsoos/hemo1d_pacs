#include <cmath>
#include <filesystem>
 
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/io/network_parser.hpp"
#include "hemo1d/physics/cfl.hpp"
#include "hemo1d/simulation.hpp"
 
using namespace hemo1d;
using Catch::Approx;
 
TEST_CASE("Simulation runs the example bifurcation network end to end", "[simulation]") {
    Network network =
        io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "simple_bifurcation.json");
 
    SimulationSettings settings;
    settings.defaultPolynomialOrder = 1;
    settings.flux = FluxKind::Hll;
    settings.useSlopeLimiter = true;
 
    Simulation sim(std::move(network), settings);
    sim.addProbe("mid_omega1", 1, 0.5);
    sim.addProbe("junction_side", 3, 1.9);
 
    CHECK(sim.time() == Approx(0.0));
 
    const Real dt = 5e-6;
    sim.run(20 * dt, dt, /*recordEvery=*/1);
 
    CHECK(sim.time() == Approx(20 * dt).margin(1e-15));
 
    for (Real a : sim.state().A) {
        REQUIRE(std::isfinite(a));
        REQUIRE(a > 0.0);
    }
    for (Real q : sim.state().Q) {
        REQUIRE(std::isfinite(q));
    }
 
    const auto& samples = sim.probeSamples("mid_omega1");
    CHECK(samples.size() == 20);
    CHECK(samples.back().time == Approx(sim.time()).margin(1e-12));
 
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "hemo1d_test_sim_probes";
    std::filesystem::create_directories(dir);
    sim.writeProbesCsv(dir.string());
    CHECK(std::filesystem::exists(dir / "mid_omega1.csv"));
    CHECK(std::filesystem::exists(dir / "junction_side.csv"));
}
 
TEST_CASE("Simulation::setUniformInitialCondition overrides the default per-vessel A0", "[simulation]") {
    Network network =
        io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "simple_bifurcation.json");
    Simulation sim(std::move(network));
 
    sim.setUniformInitialCondition(0.02, 0.2);
    for (Real a : sim.state().A) CHECK(a == Approx(0.2));
    for (Real q : sim.state().Q) CHECK(q == Approx(0.02));
}
 
TEST_CASE("Simulation::fieldSnapshot reflects every DOF's coordinate and state exactly", "[simulation]") {
    Network network =
        io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "simple_bifurcation.json");
    Simulation sim(std::move(network));
    sim.setUniformInitialCondition(0.03, 0.15);
 
    const auto snapshot = sim.fieldSnapshot();
    CHECK(snapshot.size() == static_cast<size_t>(sim.state().size()));
 
    // Cross-check against the mesh's own element/DOF layout: flattening
    // elements() in order should reproduce exactly the same (elementIndex,
    // vesselId, z) sequence and the same state values at each offset.
    size_t k = 0;
    for (const dg::Element& el : sim.mesh().elements()) {
        const std::vector<Real>& nodes = el.physicalNodes();
        const Index offset = el.dofOffset();
        for (Index i = 0; i < el.numDofs(); ++i, ++k) {
            REQUIRE(k < snapshot.size());
            CHECK(snapshot[k].elementIndex == el.flatIndex());
            CHECK(snapshot[k].vesselId == el.vesselId());
            CHECK(snapshot[k].z == Approx(nodes[i]));
            CHECK(snapshot[k].A == Approx(sim.state().A[offset + i]));
            CHECK(snapshot[k].Q == Approx(sim.state().Q[offset + i]));
        }
    }
    CHECK(k == snapshot.size());
 
    for (const auto& s : snapshot) {
        CHECK(s.A == Approx(0.15));
        CHECK(s.Q == Approx(0.03));
    }
}
 
TEST_CASE("Simulation::cflTimeStep matches physics::cflTimeStep on the live state", "[simulation]") {
    Network network =
        io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "simple_bifurcation.json");
    Simulation sim(std::move(network));
    sim.setUniformInitialCondition(0.03, 0.15);
 
    const physics::LinearElasticTubeLaw law;
    const physics::BloodFlowModel model(law, sim.network().fluid());
    const Real expected = physics::cflTimeStep(sim.state(), sim.mesh(), model, 0.8);
    CHECK(sim.cflTimeStep(0.8) == Approx(expected));
 
    // dt shrinks as the safety factor shrinks, and the helper reacts to
    // state changes rather than caching a stale value.
    CHECK(sim.cflTimeStep(0.4) == Approx(0.5 * sim.cflTimeStep(0.8)));
    const Real dtAtRest = sim.cflTimeStep(0.8);
    sim.setUniformInitialCondition(0.3, 0.15);
    CHECK(sim.cflTimeStep(0.8) != Approx(dtAtRest));
}
 
TEST_CASE("Simulation::run with dt <= 0 uses dynamic CFL-based stepping and still lands on targetTime",
          "[simulation]") {
    Network network =
        io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "simple_bifurcation.json");
    Simulation sim(std::move(network));
 
    const Real targetTime = 5e-5;
    sim.run(targetTime, /*dt=*/-1.0, /*recordEvery=*/0, /*vtkDirectory=*/"", /*vtkEvery=*/0,
            /*cflNumber=*/0.5);
 
    CHECK(sim.time() == Approx(targetTime).margin(1e-15));
    for (Real a : sim.state().A) {
        REQUIRE(std::isfinite(a));
        REQUIRE(a > 0.0);
    }
    for (Real q : sim.state().Q) REQUIRE(std::isfinite(q));
}
 