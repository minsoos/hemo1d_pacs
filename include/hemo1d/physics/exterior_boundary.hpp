#pragma once
 
#include "hemo1d/core/node.hpp"
#include "hemo1d/core/vessel.hpp"
#include "hemo1d/physics/section_state.hpp"
 
namespace hemo1d::physics{

class BloodFlowModel;

struct TerminalRow {
    Real cA = 0.0;
    Real cQ = 0.0;
    Real rhs = 0.0;
};

SectionState closeTerminal(
    VesselEnd end, const TerminalRow& physical, SectionState trace,
    SectionGradient grad, const VesselParameters& p,
    const BloodFlowModel& model, Real dt
);

SectionState solveExteriorBoundary(
    VesselEnd end, const BoundaryConditionSpec& bc, Real prescribedValue,
    SectionState trace, SectionGradient grad, const VesselParameters& p,
    const BloodFlowModel& model, Real dt
);


} // namespace hemo1d::physics