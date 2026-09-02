#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

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
    Real rho = 0.0;
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

    virtual std::vector<Real> internalState() const { return {}; }
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


class TerminalCouplingRegistry {
public:
    using Builder = std::function<std::unique_ptr<TerminalCoupling>(const std::string& paramsJson)>;

    static TerminalCouplingRegistry& instance();

    void add(std::string name, Builder builder);
    bool has(const std::string& name) const { return builders_.count(name) != 0; }

    std::unique_ptr<TerminalCoupling> build(
        const std::string& name,
        const std::string& paramsJson
    ) const;

private:
    std::unordered_map<std::string, Builder> builders_;
};

std::unique_ptr<TerminalCoupling> makeTerminalCoupling(const BoundaryConditionSpec& spec);

} // namespace hemo1d::physics