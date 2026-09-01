#include "hemo1d/physics/terminal_coupling.hpp"

#include <stdexcept>
#include <utility>

#include "hemo1d/physics/exterior_boundary.hpp"

namespace hemo1d::physics {

// ------- Prescribed Coupling ---------------

PrescribedCoupling::PrescribedCoupling(
    BoundaryConditionSpec spec, io::TimeSeries series
)
    : spec_(std::move(spec)), series_(std::move(series))
{}

SectionState PrescribedCoupling::solve(
    const TerminalInterface& iface, const BloodFlowModel& model,
    Real time, Real dt
) {
    const Real value = series_.value(time + dt);
    return solveExteriorBoundary(
        iface.end, spec_, value, iface.trace, iface.grad, iface.params, model, dt
    );
}

// ------- NonReflectingCoupling --------------

SectionState NonReflectingCoupling::solve(
    const TerminalInterface& iface, const BloodFlowModel& model, Real /*time*/, Real dt
) {
    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;
    return solveExteriorBoundary(iface.end, bc, 0.0, iface.trace, iface.grad, iface.params, model, dt);
}

// -------- CallbackCoupling -------------------

CallbackCoupling::CallbackCoupling(SolveFn SolveFn, CommitFn commitFn)
    : resolve_(std::move(SolveFn)), commit_(std::move(commitFn)) 
{
    if (!resolve_) throw std::invalid_argument("CallbackCoupling: solve callback must be set");   
}

SectionState CallbackCoupling::solve(
    const TerminalInterface& iface, const BloodFlowModel& /*model*/,
    Real time, Real dt
) {
    return resolve_(iface, time, dt);
}

void CallbackCoupling::commit(
    SectionState resolved, const TerminalInterface& iface,
    const BloodFlowModel& /*model*/, Real time, Real dt
) {
    if (commit_) commit_(resolved, iface, time, dt);
}

// -------- Factory -------------

std::unique_ptr<TerminalCoupling> makeTerminalCoupling(const BoundaryConditionSpec& spec) {
    switch (spec.type) {
        case BoundaryConditionType::NonReflecting:
            return std::make_unique<NonReflectingCoupling>();
        case BoundaryConditionType::Prescribed:
            return std::make_unique<PrescribedCoupling>(spec, io::TimeSeries::fromCsv(spec.csvFile));
    }
    throw std::invalid_argument("makeTerminalCoupling: unhandled boundary condition type");
}

} // namespace hemo1d::physics