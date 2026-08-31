#include <cmath>

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

// Cross-checks BloodFlowModel's pointwise physics against independently written
// reference formulas (papers/Master_Thesis.pdf Sect. 2.2), so a regression in
// the model's kernel is caught here rather than only through a full simulation.
TEST_CASE("BloodFlowModel matches the reference blood-flow formulas", "[blood_flow_model]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = params();

    for (Real A : {0.9 * kA0, kA0, 1.15 * kA0}) {
        for (Real Q : {-0.4, 0.0, 0.25}) {
            const SectionState u{A, Q};
            const Real vel = Q / A;
            const Real c = law.waveSpeed(A, kA0, kBeta, kRho);
            const Real cAlpha = std::sqrt(c * c + kAlpha * (kAlpha - 1.0) * vel * vel);

            // Physical flux F(U) = (Q, alpha Q^2/A + pressureFluxIntegral(A)).
            const SectionState f = model.physicalFlux(u, p);
            CHECK(f.A == Approx(Q));
            CHECK(f.Q == Approx(kAlpha * Q * Q / A + law.pressureFluxIntegral(A, kA0, kBeta, kRho)));

            // Eigenvalues lambda = alpha u +/- c_alpha.
            const Characteristics ch = model.characteristics(u, p);
            CHECK(ch.lambdaMinus == Approx(kAlpha * vel - cAlpha));
            CHECK(ch.lambdaPlus == Approx(kAlpha * vel + cAlpha));

            // Left eigenvectors: (+/- c_alpha - alpha u, 1).
            const LeftEigenvectors e = model.leftEigenvectors(u, p);
            CHECK(e.lPlus.first == Approx(cAlpha - kAlpha * vel));
            CHECK(e.lPlus.second == Approx(1.0));
            CHECK(e.lMinus.first == Approx(-cAlpha - kAlpha * vel));

            // Forward-Euler compatibility prediction with a nonzero gradient.
            const Real dAdz = 0.3, dQdz = -0.2, dt = 1e-4;
            const Real hQ = (c * c - kAlpha * vel * vel) * dAdz + 2.0 * kAlpha * vel * dQdz;
            const SectionState pred = model.compatibilityPrediction(u, {dAdz, dQdz}, p, dt);
            CHECK(pred.A == Approx(A - dt * dQdz));
            CHECK(pred.Q == Approx(Q - dt * (hQ + kKr * vel)));

            // Static and total pressure.
            CHECK(model.pressure(A, p) == Approx(law.pressure(A, kA0, kBeta)));
            CHECK(model.totalPressure(u, p) == Approx(law.pressure(A, kA0, kBeta) + 0.5 * kRho * vel * vel));
        }
    }
}
