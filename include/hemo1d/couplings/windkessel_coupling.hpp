#pragma once

#include <memory>
#include <string>
#include <vector>

#include "hemo1d/core/types.hpp"
#include "hemo1d/io/time_series.hpp"
#include "hemo1d/physics/terminal_coupling.hpp"

namespace hemo1d::couplings {

struct LumpedCompartment {
    Real resistance = 0.0;
    Real compliance = 0.0;
};

struct WindkesselParameters {
    Real R1 = -1.0;
    std::vector<LumpedCompartment> compartments;
    Real pOut = 0.0;
    std::string pOutCsv;
    Real pInit = 0.0;
    int subSteps = 1;
};


class WindkesselCoupling : public physics::TerminalCoupling {
public:
    explicit WindkesselCoupling(WindkesselParameters params);

    physics::SectionState solve(
        const physics::TerminalInterface& iface, const physics::BloodFlowModel& model,
        Real time, Real dt
    ) override;

    void commit(
        physics::SectionState resolved, const physics::TerminalInterface& iface,
        const physics::BloodFlowModel& model, Real time, Real dt
    ) override;

    const std::vector<Real>& compartmentPressures() const noexcept { return pressures_; }
    std::vector<Real> internalState() const override { return pressures_; }

    Real effectiveR1() const noexcept { return r1_; }

private:
    void ensureR1(const physics::TerminalInterface& iface, const physics::BloodFlowModel& model);
    Real pOut(Real time) const;

    WindkesselParameters params_;
    io::TimeSeries pOutSeries_;
    std::vector<Real> pressures_;
    Real r1_ = 0.0;
    bool r1Ready_ = false;
};

std::unique_ptr<physics::TerminalCoupling> makeWindkesselCoupling(WindkesselParameters params);

WindkesselParameters parseWindkesselParams(const std::string& jsonText);

} // namespace hemo1d::couplings