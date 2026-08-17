#pragma once

#include <optional>
#include <string>
#include <vector>

#include "hemo1d/core/types.hpp"
#include "hemo1d/core/vessel.hpp"

namespace hemo1d {

// This models both, the Terminal nodes, connecting with
// one other node, and the Junction, connecting with more nodes
enum class NodeKind { Terminal, Junction};

// Type of connection possible
enum class BoundaryConditionType { NonReflecting, Prescribed};

// In case it is prescribed, it is a quantity:
enum class PrescribedQuantity { FlowRate, Velocity, Area, Pressure};

struct BoundaryConditionSpec {
    // Default type and quantity
    BoundaryConditionType type = BoundaryConditionType::NonReflecting;

    PrescribedQuantity quantity = PrescribedQuantity::FlowRate;

    std::string csvFile;
};

struct VesselConnection {
    Id vesselId;
    VesselEnd end;
};

class Node{
private:
    Id id_;
    std::string name_;
    std::vector<VesselConnection> connections_;
    std::optional<BoundaryConditionSpec> boundaryCondition_;
    std::vector<Real> bifurcationAnglesRad_;

public:
    Node(Id id, std::string name, 
        std::vector<VesselConnection> connections, 
        std::optional<BoundaryConditionSpec> boundaryCondition = std::nullopt,
        std::vector<Real> bifurcationAnglesRad = {});

    Id id() const noexcept {return id_;}
    const std::string& name() const noexcept {return name_;}
    const std::vector<VesselConnection>& connections() const noexcept {return connections_;}
    const std::optional<BoundaryConditionSpec>& boundaryCondition() const noexcept {return boundaryCondition_;}

    // It has connections().size() - 1 length.
    // The angle k is between connections()[0] and connections()[k]
    const std::vector<Real>& bifurcationAnglesRad() const noexcept {return bifurcationAnglesRad_;}

    NodeKind kind() const noexcept{
        return connections_.size()<=1 ? NodeKind::Terminal : NodeKind::Junction;
    }
};

} // namespace hemo1d