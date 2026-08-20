#include <cmath>
 
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/dg/quadrature.hpp"
 
using namespace hemo1d;
using namespace hemo1d::dg;
using Catch::Approx;
 
namespace {
 
// Integrates x^power over [-1, 1] using the given rule.
double integratePower(const QuadratureRule& rule, int power) {
    double sum = 0.0;
    for (std::size_t i = 0; i < rule.points.size(); ++i) {
        sum += rule.weights[i] * std::pow(rule.points[i], power);
    }
    return sum;
}
 
// Exact value of integral_{-1}^{1} x^power dx.
double exactPowerIntegral(int power) {
    if (power % 2 != 0) return 0.0;
    return 2.0 / (power + 1);
}
 
} // namespace
 
TEST_CASE("gaussLegendre integrates polynomials exactly up to degree 2n-1", "[quadrature]") {
    for (unsigned n = 1; n <= 6; ++n) {
        const QuadratureRule rule = gaussLegendre(n);
        REQUIRE(rule.points.size() == n);
 
        Real weightSum = 0.0;
        for (Real w : rule.weights) weightSum += w;
        CHECK(weightSum == Approx(2.0));
 
        const int maxExactDegree = static_cast<int>(2 * n - 1);
        for (int power = 0; power <= maxExactDegree; ++power) {
            CHECK(integratePower(rule, power) == Approx(exactPowerIntegral(power)).margin(1e-10));
        }
    }
}
 
TEST_CASE("gaussLegendre rejects zero points", "[quadrature]") { CHECK_THROWS(gaussLegendre(0)); }
 
TEST_CASE("gaussLobattoLegendre includes the endpoints and sums weights to 2",
          "[quadrature]") {
    for (unsigned numPoints = 2; numPoints <= 7; ++numPoints) {
        const QuadratureRule rule = gaussLobattoLegendre(numPoints);
        REQUIRE(rule.points.size() == numPoints);
 
        CHECK(rule.points.front() == Approx(-1.0));
        CHECK(rule.points.back() == Approx(1.0));
 
        Real weightSum = 0.0;
        for (Real w : rule.weights) weightSum += w;
        CHECK(weightSum == Approx(2.0));
 
        const int maxExactDegree = static_cast<int>(2 * numPoints - 3);
        for (int power = 0; power <= maxExactDegree; ++power) {
            CHECK(integratePower(rule, power) == Approx(exactPowerIntegral(power)).margin(1e-10));
        }
    }
}
 
TEST_CASE("gaussLobattoLegendre matches the known 3-point rule", "[quadrature]") {
    const QuadratureRule rule = gaussLobattoLegendre(3);
    REQUIRE(rule.points.size() == 3);
    CHECK(rule.points[0] == Approx(-1.0));
    CHECK(rule.points[1] == Approx(0.0));
    CHECK(rule.points[2] == Approx(1.0));
    CHECK(rule.weights[0] == Approx(1.0 / 3.0));
    CHECK(rule.weights[1] == Approx(4.0 / 3.0));
    CHECK(rule.weights[2] == Approx(1.0 / 3.0));
}
 
TEST_CASE("gaussLobattoLegendre rejects fewer than 2 points", "[quadrature]") {
    CHECK_THROWS(gaussLobattoLegendre(1));
}