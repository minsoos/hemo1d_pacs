#pragma once

#include <utility>

#include "hemo1d/core/types.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {

// Conservative form of flux in the 1D blood flow equations.
inline std::pair<Real, Real> physicalFlux(
    Real A, Real Q, Real A0, Real beta, Real alpha, Real rho, const TubeLaw& tubeLaw
) {
    return {
        Q, 
        alpha * Q * Q / A + tubeLaw.pressureFluxIntegral(A, A0, beta, rho)
    };
}

} // namespace hemo1d::physics