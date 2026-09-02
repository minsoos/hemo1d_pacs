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

// -------- Registry ---------------

TerminalCouplingRegistry& TerminalCouplingRegistry::instance() {
    static TerminalCouplingRegistry registry;
    return registry;
}

void TerminalCouplingRegistry::add(std::string name, Builder builder) {
    if (!builder) throw std::invalid_argument(
        "TerminalCouplingRegistry::add: null builder"
    );

    auto [it, inserted] = builders_.emplace(std::move(name), std::move(builder));
    if (!inserted) {
        throw std::invalid_argument(
            "TerminalCouplingRegistry::add: model '" + it->first + 
            "' is already registered"
        );
    }
}

std::unique_ptr<TerminalCoupling> TerminalCouplingRegistry::build(
    const std::string& name, const std::string& paramsJson
) const {
    const auto it = builders_.find(name);
    if (it == builders_.end()) {
        throw std::invalid_argument(
            "makeTerminalCoupling: no terminal-coupling model named '" + name +
            "' is registered (is hemo1d_models linked and registerBuiltinModels() called?)"
        );
    }
    return it->second(paramsJson);
}

// -------- Factory -------------

std::unique_ptr<TerminalCoupling> makeTerminalCoupling(const BoundaryConditionSpec& spec) {
    switch (spec.type) {
        case BoundaryConditionType::NonReflecting:
            return std::make_unique<NonReflectingCoupling>();
        case BoundaryConditionType::Prescribed:
            return std::make_unique<PrescribedCoupling>(spec, io::TimeSeries::fromCsv(spec.csvFile));
        case BoundaryConditionType::External:
            return TerminalCouplingRegistry::instance().build(spec.modelName, spec.modelParams);
    }
    throw std::invalid_argument("makeTerminalCoupling: unhandled boundary condition type");
}

} // namespace hemo1d::physics