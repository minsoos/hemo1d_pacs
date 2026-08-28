#pragma once

#include <array>
#include <unordered_map>
#include <utility>

#include "hemo1d/core/network.hpp"
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/io/time_series.hpp"
#include "hemo1d/physics/boundary_state_provider.hpp"
#include "hemo1d/physics/state.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {
class BoundarySolver : public BoundaryStateProvider {
private:
    const Network& network_;
    const dg::Mesh& mesh_;
    const TubeLaw& tubeLaw_;
    FluidProperties fluid_;

    std::unordered_map<Id, io::TimeSeries> prescribedSeries_; // keyed by node id

    std::unordered_map<Id, std::array<std::pair<Real, Real>, 2>> ghostState_;


public:
    BoundarySolver(const Network& network, const dg::Mesh& mesh, const TubeLaw& tubeLaw,
                      FluidProperties fluid);

    void solve(const State& state, Real time, Real dt);

    std::pair<Real, Real> exteriorState(Id vesselId, VesselEnd end, Real time) const override;

};
} // namespace hemo1d::physics