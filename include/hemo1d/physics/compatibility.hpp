#pragma once

#include <utility>
 
#include "hemo1d/core/types.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {

struct LeftEigenvectors {
    std::pair<Real, Real> lPlus;
    std::pair<Real, Real> lMinus;
};

LeftEigenvectors computeLeftEigenvectors(Real A, Real Q, Real A0,
                    Real beta, Real alpha, Real rho, const TubeLaw& tubelaw);


std::pair<Real, Real> compatibilityPrediction(Real A, Real Q, Real dAdz, Real dQdz, 
                    Real A0, Real beta, Real alpha, Real rho, Real frictionKr, 
                    const TubeLaw& tubelaw, Real dt);

                    
} // namespace hemo1d::physics