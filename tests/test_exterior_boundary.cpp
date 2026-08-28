#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <tuple>

#include "hemo1d/physics/exterior_boundary.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace{

VesselParameters makeParams(){
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


TEST_CASE("solveExteriorBoundary: Prescribed: Area", "[exterior_boundary]") {
    LinearElasticTubeLaw law;
    const VesselParameters p = makeParams();

    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::Prescribed;
    bc.quantity = PrescribedQuantity::Area;
    const Real target = p.A0 * 1.02;
    const Real Q0 = 0.05;

    const auto [gA, gQ] =
        solveExteriorBoundary(VesselEnd::Proximal, bc, target, p.A0, Q0, 0.0, 0.0, p, kRho, law, kDt);
    CHECK(gA == Approx(target));
    std::ignore = gQ;
}

TEST_CASE("solveExteriorBoundary: Prescribed: Flux", "[exterior_boundary]") {
    LinearElasticTubeLaw law;
    const VesselParameters p = makeParams();

    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::Prescribed;
    bc.quantity = PrescribedQuantity::FlowRate;
    const Real target = 0.07;
    const Real Q0 = 0.05;

    const auto [gA, gQ] =
        solveExteriorBoundary(VesselEnd::Proximal, bc, target, p.A0, Q0, 0.0, 0.0, p, kRho, law, kDt);
    CHECK(gQ == Approx(target));
    std::ignore = gA;
}

TEST_CASE("solveExteriorBoundary: Prescribed: Pressure: Check compatibility areaFromPressure", "[exterior_boundary]") {
    LinearElasticTubeLaw law;
    const VesselParameters p = makeParams();

    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::Prescribed;
    bc.quantity = PrescribedQuantity::Pressure;
    const Real target = 500.0;
    const Real Q0 = 0.05;

    const auto [gA, gQ] =
        solveExteriorBoundary(VesselEnd::Proximal, bc, target, p.A0, Q0, 0.0, 0.0, p, kRho, law, kDt);
    CHECK(gA == Approx(law.areaFromPressure(target, p.A0, p.beta)));
    std::ignore = gQ;
}
