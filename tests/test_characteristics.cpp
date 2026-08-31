#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/physics/blood_flow_model.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kRho = 1.05;
constexpr Real kAlpha = 4.0 / 3.0;

const LinearElasticTubeLaw kLaw;
const BloodFlowModel kModel(kLaw, FluidProperties{kRho});

VesselParameters params() {
    VesselParameters p;
    p.A0 = kA0;
    p.beta = kBeta;
    p.alpha = kAlpha;
    return p;
}
} // namespace

TEST_CASE("characteristics at rest gives symmetric +/- wave speeds", "[characteristics]") {
    const Characteristics c = kModel.characteristics({kA0, 0.0}, params());
    const Real c0 = kLaw.waveSpeed(kA0, kA0, kBeta, kRho);
    CHECK(c.lambdaMinus == Approx(-c0));
    CHECK(c.lambdaPlus == Approx(c0));
}

TEST_CASE("characteristics orders lambdaMinus below lambdaPlus", "[characteristics]") {
    for (Real Q : {-0.5, -0.1, 0.0, 0.1, 0.5}) {
        const Characteristics c = kModel.characteristics({kA0, Q}, params());
        CHECK(c.lambdaMinus < c.lambdaPlus);
    }
}

TEST_CASE("characteristics shifts with the flow velocity", "[characteristics]") {
    const Characteristics atRest = kModel.characteristics({kA0, 0.0}, params());
    const Characteristics flowing = kModel.characteristics({kA0, 0.2}, params());
    CHECK(flowing.lambdaPlus > atRest.lambdaPlus);
    CHECK(flowing.lambdaMinus > atRest.lambdaMinus);
}
