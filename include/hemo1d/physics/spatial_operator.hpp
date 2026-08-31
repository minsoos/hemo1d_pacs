#pragma once

#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/boundary_state_provider.hpp"
#include "hemo1d/physics/flux.hpp"
#include "hemo1d/physics/state.hpp"

namespace hemo1d::physics {

class BloodFlowModel;

// Assembles the semi-discrete DG residual dU/dt = L(U) for the interior problem.
// For every element, it calculates:
//  - Volume integral of the physical flux
//  - Friction source term
//  - Numerical flux at both interfaces
class SpatialOperator {
public:
    SpatialOperator(
        const dg::Mesh& mesh, const BloodFlowModel& model, 
        const NumericalFlux& flux, const BoundaryStateProvider& boundaryProvider
    );

    void evaluate(const State& u, Real time, State& dudt) const;

    const dg::Mesh& mesh() const noexcept { return mesh_; }

private:
    const dg::Mesh& mesh_;
    const BloodFlowModel& model_;
    const NumericalFlux& flux_;
    const BoundaryStateProvider& boundaryProvider_;
};

} // namespace::hemo1d::physics