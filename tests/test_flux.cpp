#include <cmath>
#include <stdexcept>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/flux.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kRho = 1.05;
constexpr Real kAlpha = 4.0 / 3.0;

VesselParameters testParams() {
    VesselParameters p;
    p.A0 = kA0;
    p.beta = kBeta;
    p.alpha = kAlpha;
    return p;
}
} // namespace

TEMPLATE_TEST_CASE("NumericalFlux is consistent: F*(U, U) == F(U)", "[flux]", LaxFriedrichsFlux,
                    HllFlux) {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = testParams();
    TestType flux;
    for (Real Q : {-0.3, 0.0, 0.3}) {
        const SectionState f = model.physicalFlux({kA0, Q}, p);
        const SectionState star = flux.compute({kA0, Q}, {kA0, Q}, p, model);
        CHECK(star.A == Approx(f.A).margin(1e-10));
        CHECK(star.Q == Approx(f.Q).margin(1e-10));
    }
}

namespace {
// Finds a |Q| large enough that alpha*u dominates the wave speed c_alpha, so
// both characteristics at (A, Q) share the sign of Q. Avoids hard-coding a
// magic velocity that happens to work for one particular tube law.
Real dominantMomentum(Real A, const VesselParameters& p, const BloodFlowModel& model) {
    Real Q = A; // start with u = Q/A = 1
    for (int i = 0; i < 60; ++i) {
        const Characteristics c = model.characteristics({A, Q}, p);
        if (c.lambdaMinus > 0.0) return 4.0 * Q; // extra margin so both trace areas are dominated
        Q *= 2.0;
    }
    throw std::runtime_error("dominantMomentum: failed to dominate the wave speed");
}
} // namespace

TEST_CASE("HllFlux reduces to the left physical flux when both waves move right", "[flux]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = testParams();
    HllFlux flux;
    // Use the larger (harder to dominate) of the two trace areas to pick Q.
    const Real Q = dominantMomentum(kA0 * 1.05, p, model);
    const SectionState f = model.physicalFlux({kA0, Q}, p);
    const SectionState star = flux.compute({kA0, Q}, {kA0 * 1.05, Q}, p, model);
    CHECK(star.A == Approx(f.A).margin(1e-6));
    CHECK(star.Q == Approx(f.Q).margin(1e-6));
}

TEST_CASE("HllFlux reduces to the right physical flux when both waves move left", "[flux]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = testParams();
    HllFlux flux;
    const Real Q = -dominantMomentum(kA0 * 1.05, p, model);
    const SectionState f = model.physicalFlux({kA0 * 1.05, Q}, p);
    const SectionState star = flux.compute({kA0, Q}, {kA0 * 1.05, Q}, p, model);
    CHECK(star.A == Approx(f.A).margin(1e-6));
    CHECK(star.Q == Approx(f.Q).margin(1e-6));
}

TEST_CASE("LaxFriedrichsFlux adds dissipation proportional to the jump", "[flux]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = testParams();
    LaxFriedrichsFlux flux;
    // Same left/right area but a jump in Q: the flux should differ from the
    // plain average of the physical fluxes by a nonzero dissipation term.
    const SectionState avg = model.physicalFlux({kA0, 0.05}, p);
    const SectionState avg2 = model.physicalFlux({kA0, -0.05}, p);
    const SectionState star = flux.compute({kA0, 0.05}, {kA0, -0.05}, p, model);
    CHECK(std::abs(star.Q - 0.5 * (avg.Q + avg2.Q)) > 1e-6);
    (void)avg.A;
    (void)avg2.A;
    (void)star.A;
}
