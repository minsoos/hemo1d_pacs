#include "hemo1d/physics/boundary_solver.hpp"

#include <stdexcept>
#include <string>

#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/exterior_boundary.hpp"
#include "hemo1d/physics/junction_solver.hpp"
 
namespace hemo1d::physics {

namespace{

const dg::Element& boundaryElement(const dg::Mesh& mesh, Id vesselId, VesselEnd end){
    const auto [begin_, end_] = mesh.vesselElementRange(vesselId);
    return (end == VesselEnd::Proximal) ? mesh.elements()[begin_] : mesh.elements()[end_-1];
}

struct BoundaryTrace {
    SectionState u;
    SectionGradient g;
};

const BoundaryTrace traceAndDerivative(const dg::Element& el, VesselEnd end, const State& state){
    const Index n = el.numDofs(); //Dofs
    const Index offset = el.dofOffset();
    const Index nodeIdx = (end==VesselEnd::Proximal) ? 0 : n-1;
    const DenseMatrix& D = el.referenceElement().differentiationMatrix();

    Real dAdr = 0.0, dQdr = 0.0;
    for (Index j=0; j<n; ++j){
        dAdr += D(nodeIdx, j) * state.A[offset+j];
        dQdr += D(nodeIdx, j) * state.Q[offset+j];
    }
    const Real invJ = 1.0 / el.jacobian();
    return {{state.A[offset+nodeIdx], state.Q[offset+nodeIdx]}, {dAdr*invJ, dQdr*invJ}};
}

} // namespace

BoundarySolver::BoundarySolver(const Network& network, const dg::Mesh& mesh, const BloodFlowModel& model)
    : network_(network), mesh_(mesh), model_(model) {
    for (const Node& node : network_.nodes()) {
        if (node.kind() == NodeKind::Terminal &&
            node.boundaryCondition()->type == BoundaryConditionType::Prescribed) {
            prescribedSeries_.emplace(node.id(), io::TimeSeries::fromCsv(node.boundaryCondition()->csvFile));
        }
    }
}

void BoundarySolver::solve(const State& state, Real time, Real dt){
    for (const Node& node : network_.nodes()){
        if (node.kind() == NodeKind::Terminal){
            const VesselConnection& conn = node.connections().front();
            const dg::Element& el = boundaryElement(mesh_, conn.vesselId, conn.end);
            const BoundaryTrace trace = traceAndDerivative(el, conn.end, state);

            const BoundaryConditionSpec& bc = *node.boundaryCondition();

            Real prescribedValue = 0.0;
            if (bc.type == BoundaryConditionType::Prescribed){
                prescribedValue = prescribedSeries_.at(node.id()).value(time+dt);
            }

            const SectionState ghost = solveExteriorBoundary(
                conn.end, bc, prescribedValue, trace.u, trace.g,
                el.vesselParameters(), model_, dt
            );
            
            ghostState_[conn.vesselId][conn.end == VesselEnd::Proximal ? 0 : 1] = {ghost.A, ghost.Q};

        } else{
            std::vector<JunctionBranch> branches;
            branches.reserve(node.connections().size());
            for (const VesselConnection& conn: node.connections()){
                const dg::Element& el = boundaryElement(mesh_, conn.vesselId, conn.end);
                const BoundaryTrace trace = traceAndDerivative(el, conn.end, state);
                branches.push_back({el.vesselParameters(), conn.end, trace.u, trace.g});
            }

            const JunctionSolution sol = solveJunction(branches, model_, dt);
            for (Index i=0; i<node.connections().size(); ++i){
                const VesselConnection& conn = node.connections()[i];
                ghostState_[conn.vesselId][conn.end == VesselEnd::Proximal ? 0 : 1] = {sol.A[i], sol.Q[i]};
            }

        }
    }
}

std::pair<Real, Real> BoundarySolver::exteriorState(Id vesselId, VesselEnd end, Real /*time*/) const {
    const auto it = ghostState_.find(vesselId);
    if (it == ghostState_.end()){
                throw std::out_of_range("BoundaryResolver: exteriorState: vessel " + std::to_string(vesselId) +
                                 " has not been solved yet");
    }
    return it->second[end==VesselEnd::Proximal ? 0 : 1];
}


} // namespace hemo1d::physics