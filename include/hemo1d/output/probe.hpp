#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/state.hpp"

namespace hemo1d::physics {
    class BloodFlowModel;
}

namespace hemo1d::output {
 
// One recorded sample of a Probe's quantities at a point in time.
struct ProbeSample {
    Real time = 0.0;
    Real A = 0.0;
    Real Q = 0.0;
    Real pressure = 0.0;
    Real velocity = 0.0;
};

struct FieldSample {
    Index elementIndex = 0;
    Id vesselId = 0;
    Real z = 0.0;
    Real A = 0.0;
    Real Q = 0.0;
};

class Probe{
private:
    std::string name_;
    Id vesselId_;
    Real z_;
    const dg::Element* element_;
    Real referenceCoordinate_;

public:
    Probe(std::string name, Id vesselId, Real z, const dg::Mesh& mesh);
    const std::string& name() const noexcept{return name_;}
    Id vesselId() const noexcept{return vesselId_;}
    Real z() const noexcept{return z_;}

    ProbeSample sample(
        const physics::State& state, Real time, 
        const physics::BloodFlowModel& model
    ) const;
};

class ProbeRecorder{
private:
    std::vector<Probe> probes_;
    std::vector<std::vector<ProbeSample>> history_;
public:
    Probe& addProbe(std::string name, Id vesselId, Real z, const dg::Mesh& mesh);
    void record(const physics::State& state, Real time, const physics::BloodFlowModel& model);
    const std::vector<ProbeSample>& samples(const std::string& name) const;
    void writeCsv(const std::filesystem::path& directory) const;

};
 
} // namespace hemo1d::output