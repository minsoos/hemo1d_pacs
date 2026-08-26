#pragma once

#include "hemo1d/core/network.hpp"
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/state.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {

// Largest characteristic speed magnitude over every DOF; use with a target
// CFL number to pick a stable explicit time step.
Real maxWaveSpeed(
    const State& u, const dg::Mesh& mesh, 
    const FluidProperties& fluid,
    const TubeLaw& tubeLaw
);

// Largest stable explicit dt fir the current state under the standard
// DG bound
Real cflTimeStep(
    const State& u, const dg::Mesh& mesh, const FluidProperties& fluid,
    const TubeLaw& tubeLaw, Real cflNumber = 0.9
);

} // namespace hemo1d::physics
