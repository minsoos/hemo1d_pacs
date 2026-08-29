#include "hemo1d/output/probe.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>

namespace hemo1d::output {

namespace {

const dg::Element& locateElement(Id vesselId, Real z, const dg::Mesh& mesh, Real& referenceCoordinate){
    const auto [begin, end] = mesh.vesselElementRange(vesselId);
    constexpr Real kTolerance = 1e-9;
    for (Index e=begin; e<end; ++e){
        const dg::Element& el = mesh.elements()[e];
        if (z > el.zLeft() -kTolerance && el.zRight() + kTolerance){
            Real r = 2.0 * (z-el.zLeft()) / el.length() - 1.0;
            r = std::clamp(r, -1.0, 1.0);
            referenceCoordinate = r;
            return el;
        }
    }
    throw std::out_of_range("Probe: z=" + std::to_string(z) + " is outside vessel ");

}


} // namespace

Probe::Probe(std::string name, Id vesselId, Real z, const dg::Mesh& mesh):
        name_(std::move(name)), vesselId_(vesselId), z_(z),
        element_(&locateElement(vesselId, z, mesh, referenceCoordinate_)) {}

ProbeSample Probe::sample(const physics::State& state, Real time, const physics::TubeLaw& tubeLaw) const {
    const std::vector<Real> l = element_->referenceElement().basis().evaluate(referenceCoordinate_);
    const Index offset = element_->dofOffset();

    Real A = 0.0;
    Real Q = 0.0;
    for (Index i=0; i<l.size(); ++i) {
        A += l[i] * state.A[offset + i];
        Q += l[i] * state.Q[offset + i];
    }

    const VesselParameters& p = element_->vesselParameters();
    ProbeSample result;
    result.time = time;
    result.A = A;
    result.Q = Q;
    result.pressure = tubeLaw.pressure(A, p.A0, p.beta);
    result.velocity = Q / A;
    return result;
}

Probe& ProbeRecorder::addProbe(std::string name, Id vesselId, Real z, const dg::Mesh& mesh) {
    probes_.emplace_back(std::move(name), vesselId, z, mesh);
    history_.emplace_back();
    return probes_.back();
}

void ProbeRecorder::record(const physics::State& state, Real time, const physics::TubeLaw& tubeLaw) {
    for (Index i=0; i<probes_.size(); ++i) {
        history_[i].push_back(probes_[i].sample(state, time, tubeLaw));
    }
}

const std::vector<ProbeSample>& ProbeRecorder::samples(const std::string& name) const {
    for (Index i=0; i<probes_.size(); ++i){
        if (probes_[i].name() == name) return history_[i];
    }
    throw std::out_of_range("ProbeRecorder: samples: no probe named " + name);
}

void ProbeRecorder::writeCsv(const std::filesystem::path& directory) const {
    for (Index i=0; i<probes_.size(); ++i){
        std::ofstream file(directory / (probes_[i].name() + ".csv"));
        if (!file){
            throw std::runtime_error("ProbeRecorder: writeCsv: cannot open csv for '" + probes_[i].name() + "'");
        }
        file << "time,area,flowrate,pressure,velocity\n";
        for (const ProbeSample& s : history_[i]){
            file << s.time << ',' << s.A << ',' << s.Q  << ',' << s.pressure << ',' <<
                s.velocity << '\n';
        }
    }
}

} // namespace hemo1d::output