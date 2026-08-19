#pragma once 

#include "hemo1d/core/types.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {


// The two characteristic (eigenvalue) speeds of the 1D
// blood flow hyperbolic equations at a given state.
struct Characteristics {
    Real lambdaMinus;  // alpha * U - c_alpha
    Real lambdaPlus;   // alpha * U + c_alpha
};


Characteristics computeCharacteristics(
    Real A, Real Q, Real A0, Real beta, Real alpha, Real rho, const TubeLaw& tubeLaw
);

} // namespace hemo1d::physics