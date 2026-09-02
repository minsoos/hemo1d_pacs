#include "hemo1d/couplings/windkessel_coupling.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/exterior_boundary.hpp"

namespace hemo1d::couplings {

using physics::BloodFlowModel;
using physics::SectionState;
using physics::TerminalInterface;
using physics::TerminalRow;

WindkesselCoupling::WindkesselCoupling(WindkesselParameters params)
    : params_(std::move(params)), pressures_(params_.compartments.size(), params_.pInit) 
{
    if (params_.compartments.empty()) {
        throw std::invalid_argument("WindkesselCoupling: needs at least one compartment");
    }
    if (params_.R1 == 0.0) {
        throw std::invalid_argument("WindkesselCoupling: R1 == 0 (use R1 < 0 for the matched impedance)");
    }
    for (const LumpedCompartment& c : params_.compartments) {
        if (c.compliance <= 0.0 || c.resistance <= 0.0) {
            throw std::invalid_argument("WindkesselCoupling: every compartment needs R > 0 and C > 0");
        }
    }
    if (params_.subSteps < 1) {
        throw std::invalid_argument("WindkesselCoupling: subSteps must be >= 1");
    }
    if (!params_.pOutCsv.empty()) {
        pOutSeries_ = io::TimeSeries::fromCsv(params_.pOutCsv);
    }
    if (params_.R1 > 0.0) {
        r1_ = params_.R1;
        r1Ready_ = true;
    }
}

void WindkesselCoupling::ensureR1(const TerminalInterface& iface, const BloodFlowModel& model) {
    if (r1Ready_) return;

    const Real c0 = model.waveSpeed(iface.params.A0, iface.params);
    r1_ = model.density() * c0 / iface.params.A0;
    r1Ready_ = true;
}

Real WindkesselCoupling::pOut(Real time) const {
    return pOutSeries_.times().empty() ? params_.pOut : pOutSeries_.value(time);
}

SectionState WindkesselCoupling::solve(
    const TerminalInterface& iface, const BloodFlowModel& model,
    Real /*time*/, Real dt
) {
    ensureR1(iface, model);

    // Linearize the tube law about the interior trace: p(A) ~= p0 + dp*(A - A*).
    const Real aStar = iface.trace.A;
    const Real p0 = model.pressure(aStar, iface.params);
    const Real dp = model.tubeLaw().pressureDerivative(aStar, iface.params);

    const Real sgn = (iface.end == VesselEnd::Distal) ? 1.0 : -1.0;

    TerminalRow row;
    row.cA = dp;
    row.cQ = -sgn * r1_;
    row.rhs = pressures_.front() - p0 + dp * aStar;

    return physics::closeTerminal(iface.end, row, iface.trace, iface.grad, iface.params, model, dt);
}

void WindkesselCoupling::commit(
    SectionState resolved, const TerminalInterface& iface,
    const BloodFlowModel& /*model*/, Real time, Real dt
) {
    const std::size_t n = pressures_.size();
    const Real sgn = (iface.end == VesselEnd::Distal) ? 1.0 : -1.0;
    const Real qInterface = sgn * resolved.Q;
    const Real pout = pOut(time + dt);
    const int sub = params_.subSteps;
    const Real h = dt / static_cast<Real>(sub);

    std::vector<Real> qOut(n);
    for (int s = 0; s < sub; ++s) {
        for (std::size_t k = 0; k < n; ++k) {
            const Real next = (k + 1 < n) ? pressures_[k + 1] : pout;
            qOut[k] = (pressures_[k] - next) / params_.compartments[k].resistance;
        }
        for (std::size_t k = 0; k < n; ++k) {
            const Real qIn = (k == 0) ? qInterface : qOut[k - 1];
            pressures_[k] += h / params_.compartments[k].compliance * (qIn - qOut[k]);
        }
    }
}

std::unique_ptr<physics::TerminalCoupling> makeWindkesselCoupling(WindkesselParameters params) {
    return std::make_unique<WindkesselCoupling>(std::move(params));
}

WindkesselParameters parseWindkesselParams(const std::string& jsonText) {
    WindkesselParameters wk;
    if (jsonText.empty()) return wk;

    const nlohmann::json j = nlohmann::json::parse(jsonText);
    wk.R1 = j.value("r1", -1.0);
    if (j.contains("compartments")) {
        for (const nlohmann::json& c : j.at("compartments")) {
            wk.compartments.push_back(
                LumpedCompartment{c.at("r").get<Real>(), c.at("c").get<Real>()});
        }
    }
    wk.pInit = j.value("p_init", 0.0);
    wk.subSteps = j.value("sub_steps", 1);
    if (j.contains("p_out_csv")) {
        wk.pOutCsv = j.at("p_out_csv").get<std::string>();
    } else {
        wk.pOut = j.value("p_out", 0.0);
    }
    return wk;
}

} // namespace hemo1d::couplings
