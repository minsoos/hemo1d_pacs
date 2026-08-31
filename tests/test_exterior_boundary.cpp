#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/exterior_boundary.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {

VesselParameters makeParams() {
    VesselParameters p;
    p.length = 1.0;
    p.A0 = 0.126;
    p.beta = 6.06e5;
    p.alpha = 4.0 / 3.0;
    p.frictionKr = 0.0;
    return p;
}

constexpr Real kRho = 1.05;
constexpr Real kDt = 1e-4;

} // namespace

TEST_CASE("solveExteriorBoundary with a prescribed area forces that area exactly",
          "[exterior_boundary]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::Prescribed;
    bc.quantity = PrescribedQuantity::Area;
    const Real target = p.A0 * 1.02;

    const SectionState g =
        solveExteriorBoundary(VesselEnd::Proximal, bc, target, {p.A0, 0.05}, {0.0, 0.0}, p, model, kDt);
    CHECK(g.A == Approx(target));
}

TEST_CASE("solveExteriorBoundary with a prescribed flow rate forces that flow rate exactly",
          "[exterior_boundary]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::Prescribed;
    bc.quantity = PrescribedQuantity::FlowRate;
    const Real target = 0.09;

    const SectionState g =
        solveExteriorBoundary(VesselEnd::Distal, bc, target, {p.A0, 0.05}, {0.0, 0.0}, p, model, kDt);
    CHECK(g.Q == Approx(target));
}

TEST_CASE("solveExteriorBoundary with a prescribed pressure matches areaFromPressure",
          "[exterior_boundary]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::Prescribed;
    bc.quantity = PrescribedQuantity::Pressure;
    const Real targetPressure = 500.0;

    const SectionState g = solveExteriorBoundary(VesselEnd::Proximal, bc, targetPressure, {p.A0, 0.05},
                                               {0.0, 0.0}, p, model, kDt);
    CHECK(g.A == Approx(law.areaFromPressure(targetPressure, p)));
}

TEST_CASE("solveExteriorBoundary non-reflecting leaves an exact steady state unchanged",
          "[exterior_boundary]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams(); // frictionKr = 0

    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;

    for (VesselEnd end : {VesselEnd::Proximal, VesselEnd::Distal}) {
        const SectionState g =
            solveExteriorBoundary(end, bc, 0.0, {p.A0, 0.05}, {0.0, 0.0}, p, model, kDt);
        CHECK(g.A == Approx(p.A0).margin(1e-10));
        CHECK(g.Q == Approx(0.05).margin(1e-10));
    }
}
