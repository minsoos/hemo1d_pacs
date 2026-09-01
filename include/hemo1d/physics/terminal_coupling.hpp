#pragma once

#include <functional>
#include <memory>

#include "hemo1d/core/node.hpp"
#include "hemo1d/core/vessel.hpp"
#include "hemo1d/io/time_series.hpp"
#include "hemo1d/physics/section_state.hpp"

namespace hemo1d::physics {

class BloodFlowModel;

struct TerminalInterface {
    VesselEnd end = VesselEnd::Distal;
    VesselParameters params{};
    SectionState trace{};
    SectionGradient grad{};
};

class TerminalCoupling {
public:
    virtual ~TerminalCoupling() = default;

    virtual SectionState solve(
        const TerminalInterface& iface, const BloodFlowModel& model,
        Real time, Real dt
    ) = 0;

    virtual void commit(
        SectionState /*resolved*/, const TerminalInterface& /*iface*/, 
        const BloodFlowModel& /*model*/, Real /*time*/, Real /*dt*/
    ) {}
};


class PrescribedCoupling : public TerminalCoupling {
public:
    PrescribedCoupling(BoundaryConditionSpec spec, io::TimeSeries series);

    SectionState solve(
        const TerminalInterface& iface, const BloodFlowModel& model,
        Real time, Real dt
    ) override;

private:
    BoundaryConditionSpec spec_;
    io::TimeSeries series_;
};


class NonReflectingCoupling : public TerminalCoupling {
public:
    SectionState solve(
        const TerminalInterface& iface, const BloodFlowModel& model, 
        Real time, Real dt
    ) override;
};


class CallbackCoupling : public TerminalCoupling {
public:
    using SolveFn = std::function<SectionState(const TerminalInterface&, Real, Real)>;
    using CommitFn = std::function<void(SectionState, const TerminalInterface&, Real, Real)>;

    explicit CallbackCoupling(SolveFn SolveFn, CommitFn commitFn = nullptr);

    SectionState solve(
        const TerminalInterface& iface, const BloodFlowModel& model, 
        Real time, Real dt
    ) override;

    void commit(
        SectionState resolved, const TerminalInterface& iface, 
        const BloodFlowModel& model, Real time, Real dt
    ) override;

private:
    SolveFn resolve_;
    CommitFn commit_;
};

std::unique_ptr<TerminalCoupling> makeTerminalCoupling(const BoundaryConditionSpec& spec);

} // namespace hemo1d::physics