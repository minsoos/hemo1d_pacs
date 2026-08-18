#pragma once
 
#include <array>
#include <unordered_map>
#include <vector>
 
#include "hemo1d/core/node.hpp"
#include "hemo1d/core/vessel.hpp"

namespace hemo1d{

// The fluid properties are the same in every vessel
struct FluidProperties{
Real density = 1.055; // rho_f
Real viscosity = 0.045; // mu_f
};

class Network{
private:
    void buildAndValidate();

    FluidProperties fluid_;
    std::vector<Vessel> vessels_;
    std::vector<Node> nodes_;

    std::unordered_map<Id, Index> vesselIndexById_;
    std::unordered_map<Id, Index> nodeIndexById_;

    std::vector<std::array<Id, 2>> vesselEndNodeId_;

public:
    Network(FluidProperties fluid, std::vector<Vessel> vessels,
            std::vector<Node> nodes);
    
            const FluidProperties& fluid() const noexcept { return fluid_; }
            const std::vector<Vessel>& vessels() const noexcept { return vessels_; }
            const std::vector<Node>& nodes() const noexcept { return nodes_; }

            const Vessel& vessel(Id id) const;
            const Node& node(Id id) const;

            Id nodeIdAt(Id vesselId, VesselEnd end) const;

            Index vesselCount() const noexcept { return vessels_.size(); }
            Index nodeCount() const noexcept { return nodes_.size(); }
};
} // namespace hemo1d