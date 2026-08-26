#pragma once
 
#include <utility>
 
#include "hemo1d/core/node.hpp"
#include "hemo1d/core/vessel.hpp"
#include "hemo1d/physics/tube_law.hpp"
 
namespace hemo1d::physics{
    std::pair<Real, Real> solveExteriorBoundary(VesselEnd end, const BoundaryConditionSpec& bc,
                            Real prescribedValue, Real A, Real Q, Real dAdz, Real dQdz, 
                            const VesselParameters& params, Real rho, const TubeLaw& tubeLaw, Real dt);


} // namespace hemo1d::physics