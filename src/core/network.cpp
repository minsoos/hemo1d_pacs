#include "hemo1d/core/network.hpp"
 
#include <stdexcept>
#include <string>
 
namespace hemo1d {

namespace{

Index endIndex(VesselEnd end) {return end == VesselEnd::Proximal ? 0 : 1; }

const char* endName(VesselEnd end) { return end== VesselEnd::Proximal ? "proximal" : "distal"; }

} // namespace

Network::Network(FluidProperties fluid, std::vector<Vessel> vessels,
            std::vector<Node> nodes):
                fluid_(fluid), vessels_(std::move(vessels)), nodes_(std::move(nodes)) {
                    buildAndValidate();
                }

void Network::buildAndValidate(){
    vesselIndexById_.clear();
    nodeIndexById_.clear();

    for (Index i=0; i < vessels_.size(); ++i) {
        const Id id = vessels_[i].id();
        if (!vesselIndexById_.emplace(id, i).second){
            throw std::runtime_error("Network: Error vessel duplicated, id: " + std::to_string(id));
        }
    }
    for (Index i=0; i < nodes_.size(); ++i) {
        const Id id = nodes_[i].id();
        if (!nodeIndexById_.emplace(id, i).second){
            throw std::runtime_error("Network: Error node duplicated, id: " + std::to_string(id));
        }
    }

    vesselEndNodeId_.assign(vessels_.size(), {kInvalidId, kInvalidId});

    for (const Node& node : nodes_){
        if (node.connections().empty()){
            throw std::runtime_error("Network: Node " + std::to_string(node.id()) + " has no vessel connections");
        }

        const bool isTerminal = node.kind() == NodeKind::Terminal;

        if (isTerminal && !(node.boundaryCondition().has_value())){
            throw std::runtime_error("Network: Terminal node " + std::to_string(node.id()) + " has no boundary condition");
        }

        if (!isTerminal && node.boundaryCondition().has_value()){
            throw std::runtime_error("Network: Junction node " + std::to_string(node.id()) + " shouldn't have a boundary condition");
        }

        if (isTerminal && node.boundaryCondition()->type == BoundaryConditionType::Prescribed &&
            node.boundaryCondition()->csvFile.empty()){
            throw std::runtime_error("Network: Prescribed boundary condition at node " + std::to_string(node.id()) + " has no csv");
        }

        if (isTerminal && node.boundaryCondition()->type == BoundaryConditionType::External && 
            node.boundaryCondition()->modelName.empty()) {
            throw std::runtime_error(
                "Network: external boundary condition at node " + std::to_string(node.id()) + " is missing a model name"
            );
        }

        if (!node.bifurcationAnglesRad().empty() && node.bifurcationAnglesRad().size() != node.connections().size()-1){
            throw std::runtime_error("Network: Bifurcation angles at node " + std::to_string(node.id()) + " don't match the node's number");
        }

        for (const VesselConnection& conn : node.connections()){
            const auto it = vesselIndexById_.find(conn.vesselId);
            if (it == vesselIndexById_.end()){
                throw std::out_of_range("Network: Node " + std::to_string(node.id()) + " has as connection the vessel " + std::to_string(conn.vesselId) + " which doesn't exist");
            }
            // Checking if the node was already assigned
            Index endSlot = endIndex(conn.end);
            std::array<Id, 2>& ends = vesselEndNodeId_[it->second];
            if (ends[endSlot]!=kInvalidId){
                throw std::runtime_error(std::string("Network: End ") + endName(conn.end) + " of vessel " + std::to_string(conn.vesselId) + 
                    "is connected to more than 1 node");
            }
            // Assigning the ends to vesselEndNodeId created before
            ends[endIndex(conn.end)] = node.id();
        }
    }
    for (const Vessel& vessel : vessels_){
        const std::array<Id, 2>& ends = vesselEndNodeId_[vesselIndexById_.at(vessel.id())];
        if (ends[0] == kInvalidId || ends[1] == kInvalidId){
            throw std::runtime_error("Network: Vessel " + std::to_string(vessel.id()) + " has the node " + std::to_string(ends[0]) + " or " + 
            std::to_string(ends[1]) + " which it is not assigned");
        }
    }
}


const Vessel& Network::vessel(Id id) const {
    const auto it = vesselIndexById_.find(id);
    if (it == vesselIndexById_.end()){
        throw std::runtime_error("Network: Vessel " + std::to_string(id) + " is unknown "); 

    }
    return vessels_[it->second];
}

const Node& Network::node(Id id) const {
    const auto it = nodeIndexById_.find(id);
    if (it == nodeIndexById_.end()){
        throw std::runtime_error("Network: Node " + std::to_string(id) + " is not listed"); 

    }
    return nodes_[it->second];
}

Id Network::nodeIdAt(Id vesselId, VesselEnd end) const {
    const auto it = vesselIndexById_.find(vesselId);
    if (it == vesselIndexById_.end()){
        throw std::runtime_error("Network: Vessel " + std::to_string(vesselId) + " is not listed"); 
    }
    return vesselEndNodeId_[it->second][endIndex(end)];
}
} // namespace hemo1d