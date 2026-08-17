#include "hemo1d/physics/tube_law.hpp"

#include <cmath>

namespace hemo1d::physics {

Real TubeLaw::waveSpeed(Real A, Real A0, Real beta, Real rho) const {
    return std::sqrt((A / rho) * pressureDerivative(A, A0, beta));
}

Real LinearElasticTubeLaw::pressure(Real A, Real A0, Real beta) const {
    return beta * (std::sqrt(A) - std::sqrt(A0)) / A0;
}

Real LinearElasticTubeLaw::pressureDerivative(Real A, Real A0, Real beta) const {
    return beta / (2.0 * A0 * std::sqrt(A));
}

Real LinearElasticTubeLaw::pressureFluxIntegral(Real A, Real A0, Real beta, Real rho) const {
    return beta / (3.0 * rho * A0) * (std::pow(A, 1.5) - std::pow(A0, 1.5));
}

Real LinearElasticTubeLaw::areaFromPressure(Real P, Real A0, Real beta) const {
    const Real sqrtA = (P * A0 / beta) + std::sqrt(A0);
    return sqrtA * sqrtA;
}

} // namespace hemo1d::physics