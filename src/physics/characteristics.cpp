#include "hemo1d/physics/characteristics.hpp"

#include <cmath>

namespace hemo1d::physics {

Characteristics computeCharacteristics(
    Real A, Real Q, Real A0, Real beta, Real alpha, Real rho, const TubeLaw& tubeLaw
) {
    
    Real u = Q / A;
    Real c = tubeLaw.waveSpeed(A, A0, beta, rho);
    Real cAlpha = std::sqrt(c * c + alpha * (alpha - 1.0) * u * u);

    return {alpha * u - cAlpha, alpha * u + cAlpha};
}

} // namespace hemo1d::physics