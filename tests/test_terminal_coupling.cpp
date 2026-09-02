#include <filesystem>
#include <memory>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/couplings/register.hpp"
#include "hemo1d/couplings/windkessel_coupling.hpp"
#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/exterior_boundary.hpp"
#include "hemo1d/physics/terminal_coupling.hpp"
#include "hemo1d/physics/tube_law.hpp"

using namespace hemo1d;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {

constexpr Real kRho = 1.05;
constexpr Real kDt = 1e-4;

VesselParameters makeParams() {
    VesselParameters p;
    p.length = 1.0;
    p.A0 = 0.126;
    p.beta = 6.06e5;
    p.alpha = 4.0 / 3.0;
    return p;
}

TerminalInterface makeInterface(const VesselParameters& p, VesselEnd end) {
    TerminalInterface iface;
    iface.end = end;
    iface.params = p;
    iface.trace = {p.A0, 0.05};
    iface.grad = {0.01, -0.02};
    return iface;
}

} // namespace

TEST_CASE("NonReflectingCoupling matches solveExteriorBoundary exactly", "[terminal_coupling]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();

    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;

    for (VesselEnd end : {VesselEnd::Proximal, VesselEnd::Distal}) {
        const TerminalInterface iface = makeInterface(p, end);
        NonReflectingCoupling coupling;
        const SectionState viaCoupling = coupling.solve(iface, model, 0.0, kDt);
        const SectionState viaFree =
            solveExteriorBoundary(end, bc, 0.0, iface.trace, iface.grad, p, model, kDt);
        CHECK(viaCoupling.A == Approx(viaFree.A));
        CHECK(viaCoupling.Q == Approx(viaFree.Q));
    }
}

TEST_CASE("makeTerminalCoupling builds the right kind and rejects unknowns", "[terminal_coupling]") {
    BoundaryConditionSpec nr;
    nr.type = BoundaryConditionType::NonReflecting;
    CHECK(dynamic_cast<NonReflectingCoupling*>(makeTerminalCoupling(nr).get()) != nullptr);

    BoundaryConditionSpec pr;
    pr.type = BoundaryConditionType::Prescribed;
    pr.quantity = PrescribedQuantity::FlowRate;
    pr.csvFile = (std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "data" / "inlet_flow_sine.csv").string();
    CHECK(dynamic_cast<PrescribedCoupling*>(makeTerminalCoupling(pr).get()) != nullptr);

    // External: dispatched through the registry by model name.
    BoundaryConditionSpec unknown;
    unknown.type = BoundaryConditionType::External;
    unknown.modelName = "no_such_model";
    CHECK_THROWS(makeTerminalCoupling(unknown));

    hemo1d::couplings::registerBuiltinCouplings();
    BoundaryConditionSpec wk;
    wk.type = BoundaryConditionType::External;
    wk.modelName = "windkessel";
    wk.modelParams = R"({"compartments":[{"r":5e4,"c":8e-7}]})";
    CHECK(dynamic_cast<hemo1d::couplings::WindkesselCoupling*>(makeTerminalCoupling(wk).get()) != nullptr);
}

TEST_CASE("CallbackCoupling forwards to the user function", "[terminal_coupling]") {
    const LinearElasticTubeLaw law;
    const BloodFlowModel model(law, FluidProperties{kRho});
    const VesselParameters p = makeParams();
    const TerminalInterface iface = makeInterface(p, VesselEnd::Distal);

    int calls = 0;
    CallbackCoupling coupling(
        [&](const TerminalInterface& in, Real time, Real dt) -> SectionState {
            ++calls;
            CHECK(time == Approx(0.25));
            CHECK(dt == Approx(kDt));
            return {in.trace.A * 1.1, 0.0}; // "rigid wall at 1.1 A*"
        });

    const SectionState g = coupling.solve(iface, model, 0.25, kDt);
    CHECK(calls == 1);
    CHECK(g.A == Approx(p.A0 * 1.1));
    CHECK(g.Q == Approx(0.0));

    // A null solve callback is rejected at construction.
    CHECK_THROWS(CallbackCoupling(nullptr));
}