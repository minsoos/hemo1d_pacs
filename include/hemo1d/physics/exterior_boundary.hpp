#pragma once
 
#include "hemo1d/core/node.hpp"
#include "hemo1d/core/vessel.hpp"
#include "hemo1d/physics/section_state.hpp"
 
namespace hemo1d::physics{

class BloodFlowModel;

SectionState solveExteriorBoundary(
    VesselEnd end, const BoundaryConditionSpec& bc, Real prescribedValue,
    SectionState trace, SectionGradient grad, const VesselParameters& p,
    const BloodFlowModel& model, Real dt
);


} // namespace hemo1d::physics