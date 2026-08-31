#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <utility>

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

VesselParameters params(Real frictionKr = 0.0) {
    VesselParameters p;
    p.A0 = kA0;
    p.beta = kBeta;
    p.alpha = kAlpha;
    p.frictionKr = frictionKr;
    return p;
}
} // namespace

TEST_CASE("leftEigenvectors satisfies l^T H = lambda l^T", "[compatibility]") {
    const Real A = kA0 * 1.1;
    const Real Q = 0.08;
    const Real u = Q / A;
    const Real c = kLaw.waveSpeed(A, kA0, kBeta, kRho);

    // H(U) = [[0, 1], [c^2 - alpha*u^2, 2*alpha*u]] (papers/Master_Thesis.pdf).
    const Real H00 = 0.0, H01 = 1.0;
    const Real H10 = c * c - kAlpha * u * u, H11 = 2.0 * kAlpha * u;

    const LeftEigenvectors eig = kModel.leftEigenvectors({A, Q}, params());
    const Characteristics ch = kModel.characteristics({A, Q}, params());

    for (const auto& [l, lambda] : 
         {std::make_pair(eig.lPlus, ch.lambdaPlus), std::make_pair(eig.lMinus, ch.lambdaMinus)}) {
        const Real lhs0 = l.first * H00 + l.second * H10;
        const Real lhs1 = l.first * H01 + l.second * H11;
        CHECK(lhs0 == Approx(lambda * l.first).margin(1e-8));
        CHECK(lhs1 == Approx(lambda * l.second).margin(1e-8));
    }
}

TEST_CASE("compatibilityPrediction is a no-op with zero derivatives and no friction",
          "[compatibility]") {
    const SectionState pred = kModel.compatibilityPrediction({kA0, 0.05}, {0.0, 0.0}, params(), 1e-4);
    CHECK(pred.A == Approx(kA0));
    CHECK(pred.Q == Approx(0.05));
}

TEST_CASE("compatibilityPrediction applies friction as expected", "[compatibility]") {
    const Real Q = 0.05;
    const Real Kr = 3.0;
    const Real dt = 1e-4;
    const SectionState pred = kModel.compatibilityPrediction({kA0, Q}, {0.0, 0.0}, params(Kr), dt);
    CHECK(pred.A == Approx(kA0));
    CHECK(pred.Q == Approx(Q - dt * Kr * Q / kA0));
}
