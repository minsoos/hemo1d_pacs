#include "hemo1d/core/node.hpp"

#include <utility>

namespace hemo1d {

Node::Node(Id id, std::string name, 
        std::vector<VesselConnection> connections, 
        std::optional<BoundaryConditionSpec> boundaryCondition,
        std::vector<Real> bifurcationAnglesRad):
            id_(id),
            name_(std::move(name)),
            connections_(std::move(connections)),
            boundaryCondition_(std::move(boundaryCondition)),
            bifurcationAnglesRad_(std::move(bifurcationAnglesRad)){}
} //namespace hemo1d