#include <cmath>
#include <stdexcept>
 
#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/physics/characteristics.hpp"
#include "hemo1d/physics/conservation_law.hpp"
#include "hemo1d/physics/flux.hpp"
 
using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;
 
namespace {
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kRho = 1.05;
constexpr Real kAlpha = 4.0 / 3.0;
} // namespace
 
TEMPLATE_TEST_CASE("NumericalFlux is consistent: F*(U, U) == F(U)", "[flux]", LaxFriedrichsFlux,
                    HllFlux) {
    LinearElasticTubeLaw law;
    TestType flux;
    for (Real Q : {-0.3, 0.0, 0.3}) {
        const auto [fA, fQ] = physicalFlux(kA0, Q, kA0, kBeta, kAlpha, kRho, law);
        const auto [starA, starQ] = flux.compute(kA0, Q, kA0, Q, kA0, kBeta, kAlpha, kRho, law);
        CHECK(starA == Approx(fA).margin(1e-10));
        CHECK(starQ == Approx(fQ).margin(1e-10));
    }
}
 
namespace {
// Finds a |Q| large enough that alpha*u dominates the wave speed c_alpha, so
// both characteristics at (A, Q) share the sign of Q. Avoids hard-coding a
// magic velocity that happens to work for one particular tube law.
Real dominantMomentum(Real A, Real A0, Real beta, Real alpha, Real rho, const TubeLaw& law) {
    Real Q = A; // start with u = Q/A = 1
    for (int i = 0; i < 60; ++i) {
        const Characteristics c = computeCharacteristics(A, Q, A0, beta, alpha, rho, law);
        if (c.lambdaMinus > 0.0) return 4.0 * Q; // extra margin so both trace areas are dominated
        Q *= 2.0;
    }
    throw std::runtime_error("dominantMomentum: failed to dominate the wave speed");
}
} // namespace
 
TEST_CASE("HllFlux reduces to the left physical flux when both waves move right", "[flux]") {
    LinearElasticTubeLaw law;
    HllFlux flux;
    // Use the larger (harder to dominate) of the two trace areas to pick Q.
    const Real Q = dominantMomentum(kA0 * 1.05, kA0, kBeta, kAlpha, kRho, law);
    const auto [fA, fQ] = physicalFlux(kA0, Q, kA0, kBeta, kAlpha, kRho, law);
    const auto [starA, starQ] = flux.compute(kA0, Q, kA0 * 1.05, Q, kA0, kBeta, kAlpha, kRho, law);
    CHECK(starA == Approx(fA).margin(1e-6));
    CHECK(starQ == Approx(fQ).margin(1e-6));
}
 
TEST_CASE("HllFlux reduces to the right physical flux when both waves move left", "[flux]") {
    LinearElasticTubeLaw law;
    HllFlux flux;
    const Real Q = -dominantMomentum(kA0 * 1.05, kA0, kBeta, kAlpha, kRho, law);
    const auto [fA, fQ] = physicalFlux(kA0 * 1.05, Q, kA0, kBeta, kAlpha, kRho, law);
    const auto [starA, starQ] = flux.compute(kA0, Q, kA0 * 1.05, Q, kA0, kBeta, kAlpha, kRho, law);
    CHECK(starA == Approx(fA).margin(1e-6));
    CHECK(starQ == Approx(fQ).margin(1e-6));
}
 
TEST_CASE("LaxFriedrichsFlux adds dissipation proportional to the jump", "[flux]") {
    LinearElasticTubeLaw law;
    LaxFriedrichsFlux flux;
    // Same left/right area but a jump in Q: the flux should differ from the
    // plain average of the physical fluxes by a nonzero dissipation term.
    const auto [avgA, avgQ] = physicalFlux(kA0, 0.05, kA0, kBeta, kAlpha, kRho, law);
    const auto [avgA2, avgQ2] = physicalFlux(kA0, -0.05, kA0, kBeta, kAlpha, kRho, law);
    const auto [starA, starQ] = flux.compute(kA0, 0.05, kA0, -0.05, kA0, kBeta, kAlpha, kRho, law);
    CHECK(std::abs(starQ - 0.5 * (avgQ + avgQ2)) > 1e-6);
    (void)avgA;
    (void)avgA2;
    (void)starA;
}