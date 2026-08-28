#include "hemo1d/physics/boundary_solver.hpp"

#include <stdexcept>
#include <string>
 
#include "hemo1d/physics/exterior_boundary.hpp"
#include "hemo1d/physics/junction_solver.hpp"
 
namespace hemo1d::physics {

BoundarySolver::BoundarySolver(const Network& network, const dg::Mesh& mesh,
                const TubeLaw& tubeLaw, FluidProperties fluid)
    : network_(network), mesh_(mesh), tubeLaw_(tubeLaw), fluid_(fluid) {
    for (const Node& node : network_.nodes()) {
        if (node.kind() == NodeKind::Terminal &&
            node.boundaryCondition()->type == BoundaryConditionType::Prescribed) {
            prescribedSeries_.emplace(node.id(), io::TimeSeries::fromCsv(node.boundaryCondition()->csvFile));
        }
    }
}

} // namespace hemo1d::physics