#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/core/types.hpp"
#include "hemo1d/physics/junction_solver.hpp"


using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace{
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kRho = 1.05;
constexpr Real kAlpha = 4.0 / 3.0;
constexpr Real kDt = 1e-4;

JunctionBranch makeBranch(VesselEnd end, Real A, Real Q) {
    JunctionBranch b;
    b.A0 = kA0;
    b.beta = kBeta;
    b.alpha = kAlpha;
    b.frictionKr = 0.0;
    b.end = end;
    b.A = A;
    b.Q = Q;
    b.dAdz = 0.0;
    b.dQdz = 0.0;
    return b;
}
} // namespace

TEST_CASE("solveJunction: Simplest case", "[junction_solver]") {
    const Real Q = 0.05;
    std::vector<JunctionBranch> branches{
        makeBranch(VesselEnd::Distal, kA0, Q),
        makeBranch(VesselEnd::Proximal, kA0, Q),
    };

    const JunctionSolution sol = solveJunction(branches, kRho, LinearElasticTubeLaw{}, kDt);

    CHECK(sol.iterations < 50);
    CHECK(sol.residualNorm < 1e-6);
    CHECK(sol.A[0] == Approx(kA0).margin(1e-9));
    CHECK(sol.A[1] == Approx(kA0).margin(1e-9));
    CHECK(sol.Q[0] == Approx(Q).margin(1e-9));
    CHECK(sol.Q[1] == Approx(Q).margin(1e-9));
}

TEST_CASE("solveJunction: mass conservation and pressure continuity", "[junction_solver]") {
    LinearElasticTubeLaw law;
    std::vector<JunctionBranch> branches{
        makeBranch(VesselEnd::Distal, kA0 * 1.01, 0.10),
        makeBranch(VesselEnd::Proximal, kA0, 0.04),
        makeBranch(VesselEnd::Proximal, kA0 * 0.99, 0.055),
    };

    const JunctionSolution sol = solveJunction(branches, kRho, law, kDt);
    CHECK(sol.iterations < 50);
    CHECK(sol.residualNorm < 1e-3);

    const Real massResidual = sol.Q[0] - sol.Q[1] - sol.Q[2];
    CHECK(massResidual == Approx(0.0).margin(1e-6));

    auto totalPressure = [&](std::size_t i) {
        const Real u = sol.Q[i] / sol.A[i];
        return law.pressure(sol.A[i], branches[i].A0, branches[i].beta) + 0.5 * kRho * u * u;
    };
    const Real p0 = totalPressure(0);
    CHECK(totalPressure(1) == Approx(p0).margin(1e-6));
    CHECK(totalPressure(2) == Approx(p0).margin(1e-6));
}

TEST_CASE("solveJunction: Rejects less than two branches", "[junction_solver]") {
    std::vector<JunctionBranch> branches{makeBranch(VesselEnd::Distal, kA0, 0.05)};
    CHECK_THROWS(solveJunction(branches, kRho, LinearElasticTubeLaw{}, kDt));
}