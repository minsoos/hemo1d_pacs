#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/junction_solver.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kRho = 1.05;
constexpr Real kAlpha = 4.0 / 3.0;
constexpr Real kDt = 1e-4;

const LinearElasticTubeLaw kLaw;
const BloodFlowModel kModel(kLaw, FluidProperties{kRho});

JunctionBranch makeBranch(VesselEnd end, Real A, Real Q) {
    JunctionBranch b;
    b.params.A0 = kA0;
    b.params.beta = kBeta;
    b.params.alpha = kAlpha;
    b.params.frictionKr = 0.0;
    b.end = end;
    b.trace = {A, Q};
    b.grad = {0.0, 0.0};
    return b;
}
} // namespace

TEST_CASE("solveJunction reproduces a steady straight-pipe fixed point", "[junction_solver]") {
    // Two identical vessels, one feeding in (distal end here) and one
    // continuing (proximal end here), both already at the same steady state.
    const Real Q = 0.05;
    std::vector<JunctionBranch> branches{
        makeBranch(VesselEnd::Distal, kA0, Q),
        makeBranch(VesselEnd::Proximal, kA0, Q),
    };

    const JunctionSolution sol = solveJunction(branches, kModel, kDt);

    CHECK(sol.iterations < 50);
    CHECK(sol.residualNorm < 1e-6);
    CHECK(sol.A[0] == Approx(kA0).margin(1e-9));
    CHECK(sol.A[1] == Approx(kA0).margin(1e-9));
    CHECK(sol.Q[0] == Approx(Q).margin(1e-9));
    CHECK(sol.Q[1] == Approx(Q).margin(1e-9));
}

TEST_CASE("solveJunction satisfies mass conservation and pressure continuity for a bifurcation",
          "[junction_solver]") {
    // One inflow (distal end), two outflows (proximal ends), a mild
    // perturbation away from the trivial steady split.
    std::vector<JunctionBranch> branches{
        makeBranch(VesselEnd::Distal, kA0 * 1.01, 0.10),
        makeBranch(VesselEnd::Proximal, kA0, 0.04),
        makeBranch(VesselEnd::Proximal, kA0 * 0.99, 0.055),
    };

    const JunctionSolution sol = solveJunction(branches, kModel, kDt);
    REQUIRE(sol.iterations < 50);
    // The stopping criterion is a relative increment, not the residual norm
    // directly, so only require the residual to have collapsed by orders of
    // magnitude from a typical O(1) initial residual -- the physically
    // meaningful checks below (mass, pressure) use a tight margin instead.
    CHECK(sol.residualNorm < 1e-3);

    // Mass conservation: inflow (distal) equals the sum of outflows (proximal).
    const Real massResidual = sol.Q[0] - sol.Q[1] - sol.Q[2];
    CHECK(massResidual == Approx(0.0).margin(1e-6));

    // Continuity of total pressure across all branches.
    auto totalPressure = [&](std::size_t i) {
        const Real u = sol.Q[i] / sol.A[i];
        return kLaw.pressure(sol.A[i], branches[i].params.A0, branches[i].params.beta) + 0.5 * kRho * u * u;
    };
    const Real p0 = totalPressure(0);
    CHECK(totalPressure(1) == Approx(p0).margin(1e-6));
    CHECK(totalPressure(2) == Approx(p0).margin(1e-6));
}

TEST_CASE("solveJunction rejects fewer than two branches", "[junction_solver]") {
    std::vector<JunctionBranch> branches{makeBranch(VesselEnd::Distal, kA0, 0.05)};
    CHECK_THROWS(solveJunction(branches, kModel, kDt));
}
