#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/physics/tube_law.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kRho = 1.05;

// A tube law reads only the wall fields (A0, beta) of VesselParameters.
VesselParameters wall(Real A0 = kA0, Real beta = kBeta) {
    VesselParameters p;
    p.A0 = A0;
    p.beta = beta;
    return p;
}
} // namespace

TEST_CASE("LinearElasticTubeLaw gives zero pressure at the reference area", "[tube_law]") {
    LinearElasticTubeLaw law;
    CHECK(law.pressure(kA0, wall()) == Approx(0.0).margin(1e-12));
    CHECK(law.pressureFluxIntegral(kA0, wall(), kRho) == Approx(0.0).margin(1e-12));
}

TEST_CASE("LinearElasticTubeLaw pressure increases with area", "[tube_law]") {
    LinearElasticTubeLaw law;
    CHECK(law.pressure(kA0 * 1.1, wall()) > 0.0);
    CHECK(law.pressure(kA0 * 0.9, wall()) < 0.0);
}

TEST_CASE("LinearElasticTubeLaw areaFromPressure inverts pressure", "[tube_law]") {
    LinearElasticTubeLaw law;
    for (Real P : {-500.0, 0.0, 250.0, 4000.0}) {
        const Real A = law.areaFromPressure(P, wall());
        CHECK(law.pressure(A, wall()) == Approx(P).margin(1e-8));
    }
}

TEST_CASE("LinearElasticTubeLaw wave speed matches the closed-form reference value",
          "[tube_law]") {
    LinearElasticTubeLaw law;
    // c(A0) = sqrt(beta / (2 rho sqrt(A0))), from c^2 = beta*sqrt(A)/(2 rho A0).
    const Real expected = std::sqrt(kBeta / (2.0 * kRho * std::sqrt(kA0)));
    CHECK(law.waveSpeed(kA0, wall(), kRho) == Approx(expected));
}

TEST_CASE("LinearElasticTubeLaw pressureFluxIntegral derivative matches c^2", "[tube_law]") {
    // d/dA pressureFluxIntegral(A) should equal c(A)^2 = (A/rho) pressureDerivative(A).
    LinearElasticTubeLaw law;
    const Real A = kA0 * 1.2;
    const Real h = 1e-6;
    const Real numericalDerivative =
        (law.pressureFluxIntegral(A + h, wall(), kRho) - law.pressureFluxIntegral(A - h, wall(), kRho)) /
        (2.0 * h);
    const Real cSquared = (A / kRho) * law.pressureDerivative(A, wall());
    CHECK(numericalDerivative == Approx(cSquared).margin(1e-6));
}
