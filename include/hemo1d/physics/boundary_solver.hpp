#pragma once

#include <array>
#include <unordered_map>
#include <utility>

#include "hemo1d/core/network.hpp"
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/io/time_series.hpp"
#include "hemo1d/physics/boundary_state_provider.hpp"
#include "hemo1d/physics/state.hpp"

namespace hemo1d::physics {

class BloodFlowModel;

class BoundarySolver : public BoundaryStateProvider {
private:
    const Network& network_;
    const dg::Mesh& mesh_;
    const BloodFlowModel& model_;

    std::unordered_map<Id, io::TimeSeries> prescribedSeries_; // keyed by node id

    std::unordered_map<Id, std::array<std::pair<Real, Real>, 2>> ghostState_;


public:
    BoundarySolver(const Network& network, const dg::Mesh& mesh, const BloodFlowModel& model);

    void solve(const State& state, Real time, Real dt);

    std::pair<Real, Real> exteriorState(Id vesselId, VesselEnd end, Real time) const override;

};
} // namespace hemo1d::physics