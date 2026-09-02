#pragma once

#include <array>
#include <unordered_map>
#include <utility>
#include <memory>
#include <vector>

#include "hemo1d/core/network.hpp"
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/boundary_state_provider.hpp"
#include "hemo1d/physics/state.hpp"
#include "hemo1d/physics/terminal_coupling.hpp"

namespace hemo1d::physics {

class BloodFlowModel;

class BoundarySolver : public BoundaryStateProvider {
private:
    const Network& network_;
    const dg::Mesh& mesh_;
    const BloodFlowModel& model_;

    std::unordered_map<Id, std::unique_ptr<TerminalCoupling>> terminalCouplings_;

    struct TerminalRecord {
        Id nodeId = kInvalidId;
        TerminalInterface iface;
        SectionState resolved;
    };
    std::vector<TerminalRecord> terminalRecords_;

    std::unordered_map<Id, std::array<std::pair<Real, Real>, 2>> ghostState_;

public:
    BoundarySolver(const Network& network, const dg::Mesh& mesh, const BloodFlowModel& model);

    void solve(const State& state, Real time, Real dt);

    void commit(const State& state, Real time, Real dt);

    void setCoupling(Id nodeId, std::unique_ptr<TerminalCoupling> coupling);

    const TerminalCoupling* coupling(Id nodeId) const;

    std::pair<Real, Real> exteriorState(Id vesselId, VesselEnd end, Real time) const override;

};
} // namespace hemo1d::physics