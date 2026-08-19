#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/physics/characteristics.hpp"
 
using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;
 
namespace {
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kRho = 1.05;
constexpr Real kAlpha = 4.0 / 3.0;
} // namespace
 
TEST_CASE("computeCharacteristics at rest gives symmetric +/- wave speeds", "[characteristics]") {
    LinearElasticTubeLaw law;
    const Characteristics c = computeCharacteristics(kA0, 0.0, kA0, kBeta, kAlpha, kRho, law);
    const Real c0 = law.waveSpeed(kA0, kA0, kBeta, kRho);
    CHECK(c.lambdaMinus == Approx(-c0));
    CHECK(c.lambdaPlus == Approx(c0));
}
 
TEST_CASE("computeCharacteristics orders lambdaMinus below lambdaPlus", "[characteristics]") {
    LinearElasticTubeLaw law;
    for (Real Q : {-0.5, -0.1, 0.0, 0.1, 0.5}) {
        const Characteristics c = computeCharacteristics(kA0, Q, kA0, kBeta, kAlpha, kRho, law);
        CHECK(c.lambdaMinus < c.lambdaPlus);
    }
}
 
TEST_CASE("computeCharacteristics shifts with the flow velocity", "[characteristics]") {
    LinearElasticTubeLaw law;
    const Characteristics atRest = computeCharacteristics(kA0, 0.0, kA0, kBeta, kAlpha, kRho, law);
    const Characteristics flowing = computeCharacteristics(kA0, 0.2, kA0, kBeta, kAlpha, kRho, law);
    CHECK(flowing.lambdaPlus > atRest.lambdaPlus);
    CHECK(flowing.lambdaMinus > atRest.lambdaMinus);
}