#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/characteristics.hpp"
#include "hemo1d/physics/compatibility.hpp"
#include "hemo1d/physics/conservation_law.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kRho = 1.05;
constexpr Real kAlpha = 4.0 / 3.0;
constexpr Real kKr = 22.0;

VesselParameters params() {
    VesselParameters p;
    p.A0 = kA0;
    p.beta = kBeta;
    p.alpha = kAlpha;
    p.frictionKr = kKr;
    return p;
}
} // namespace

// BloodFlowModel is introduced as a thin wrapper over the existing free
// functions; this pins that it reproduces them exactly, before any caller is
// migrated onto it.
TEST_CASE("BloodFlowModel reproduces the free-function physics", "[blood_flow_model]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = params();

    for (Real A : {0.9 * kA0, kA0, 1.15 * kA0}) {
        for (Real Q : {-0.4, 0.0, 0.25}) {
            const SectionState u{A, Q};

            const auto [fA, fQ] = physicalFlux(A, Q, kA0, kBeta, kAlpha, kRho, law);
            const SectionState f = model.physicalFlux(u, p);
            CHECK(f.A == Approx(fA));
            CHECK(f.Q == Approx(fQ));

            const Characteristics cRef = computeCharacteristics(A, Q, kA0, kBeta, kAlpha, kRho, law);
            const Characteristics c = model.characteristics(u, p);
            CHECK(c.lambdaMinus == Approx(cRef.lambdaMinus));
            CHECK(c.lambdaPlus == Approx(cRef.lambdaPlus));

            const LeftEigenvectors eRef = computeLeftEigenvectors(A, Q, kA0, kBeta, kAlpha, kRho, law);
            const LeftEigenvectors e = model.leftEigenvectors(u, p);
            CHECK(e.lPlus.first == Approx(eRef.lPlus.first));
            CHECK(e.lMinus.first == Approx(eRef.lMinus.first));

            const auto [pA, pQ] =
                compatibilityPrediction(A, Q, 0.3, -0.2, kA0, kBeta, kAlpha, kRho, kKr, law, 1e-4);
            const SectionState pred = model.compatibilityPrediction(u, {0.3, -0.2}, p, 1e-4);
            CHECK(pred.A == Approx(pA));
            CHECK(pred.Q == Approx(pQ));

            CHECK(model.pressure(A, p) == Approx(law.pressure(A, kA0, kBeta)));
            CHECK(model.totalPressure(u, p) ==
                  Approx(law.pressure(A, kA0, kBeta) + 0.5 * kRho * (Q / A) * (Q / A)));
        }
    }
}