#include <cmath>
#include <filesystem>
#include <string>
#include <iostream>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/io/network_parser.hpp"
#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/boundary_solver.hpp"
#include "hemo1d/physics/solver.hpp"

using namespace hemo1d;
using namespace hemo1d::dg;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {
 
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;

Network makeSingleVesselNetwork() {
    VesselParameters params;
    params.length = 1.0;
    params.A0 = kA0;
    params.beta = kBeta;
    params.alpha = 4.0 / 3.0;
    params.frictionKr = 0.0;
    params.nElements = 5;
    params.polynomialOrder = 2;

    std::vector<Vessel> vessels{Vessel(1, "v", params)};

    BoundaryConditionSpec inlet;
    inlet.type = BoundaryConditionType::Prescribed;
    inlet.quantity = PrescribedQuantity::Area;
    inlet.csvFile = (std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "data" / "inlet_area_ramp.csv").string();

    BoundaryConditionSpec outlet;
    outlet.type = BoundaryConditionType::NonReflecting;

    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, inlet),
        Node(2, "out", {{1, VesselEnd::Distal}}, outlet),
    };
    return Network(FluidProperties{}, std::move(vessels), std::move(nodes));
}

bool allFinite(const State& s) {
    for (Real v : s.A) {
        if (!std::isfinite(v) || v <= 0.0) return false;
    }
    for (Real v : s.Q) {
        if (!std::isfinite(v)) return false;
    }
    return true;
}

} // namespace

TEST_CASE("BoundarySolver: Single vessel with a prescribed inlet and stays stable",
          "[boundary_solver]") {
    const Network network = makeSingleVesselNetwork();
    const Mesh mesh(network);
 
    State state(mesh.totalDofs());
    for (Index i = 0; i < state.size(); ++i) {
        state.A[i] = kA0;
        state.Q[i] = 0.0;
    }
 
    LinearElasticTubeLaw law;
    LaxFriedrichsFlux flux;
    MinmodLimiter limiter;
    const BloodFlowModel model(law, FluidProperties{});
    BoundarySolver bSolver(network, mesh, model);
    Solver solver(mesh, model, flux, bSolver, &limiter);

    Real time = 0.0;
    const Real dt = 1e-5;
    for (int step = 0; step < 30; ++step) {
        bSolver.solve(state, time, dt);
        solver.step(state, time, dt);
        time += dt;
        REQUIRE(allFinite(state));
    }

    for (Real A : state.A) {
        CHECK(A == Approx(kA0).epsilon(0.05));
    }
}

TEST_CASE("BoundarySolver: Example implementtion",
          "[boundary_solver]") {
    const Network network =
        io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "simple_bifurcation.json");
    DgSettings settings;
    settings.defaultPolynomialOrder = 1;
    const Mesh mesh(network, settings);
    
    State state(mesh.totalDofs());
    for (const Element& el : mesh.elements()) {
        for (Index i = 0; i < el.numDofs(); ++i) {
            state.A[el.dofOffset() + i] = el.vesselParameters().A0;
            state.Q[el.dofOffset() + i] = 0.0;
        }
    }

    LinearElasticTubeLaw law;
    HllFlux flux;
    MinmodLimiter limiter;
    const BloodFlowModel model(law, FluidProperties{});
    BoundarySolver bSolver(network, mesh, model);
    Solver solver(mesh, model, flux, bSolver, &limiter);

    Real time = 0.0;
    const Real dt = 5e-6;
    for (int step = 0; step < 20; ++step) {
        bSolver.solve(state, time, dt);
        solver.step(state, time, dt);
        time += dt;
        REQUIRE(allFinite(state));
    }
}