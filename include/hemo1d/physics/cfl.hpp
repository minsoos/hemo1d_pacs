#pragma once

#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/state.hpp"

namespace hemo1d::physics {

class BloodFlowModel;

// Largest characteristic speed magnitude over every DOF; use with a target
// CFL number to pick a stable explicit time step.
Real maxWaveSpeed(
    const State& u, const dg::Mesh& mesh, const BloodFlowModel& model
);

// Largest stable explicit dt fir the current state under the standard
// DG bound
Real cflTimeStep(
    const State& u, const dg::Mesh& mesh, const BloodFlowModel& model,
    Real cflNumber = 0.9
);

} // namespace hemo1d::physics
