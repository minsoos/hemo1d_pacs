#include "hemo1d/physics/tube_law.hpp"

#include <cmath>

namespace hemo1d::physics {

Real TubeLaw::waveSpeed(Real A, const VesselParameters& p, Real rho) const {
    return std::sqrt((A / rho) * pressureDerivative(A, p));
}

Real LinearElasticTubeLaw::pressure(Real A, const VesselParameters& p) const {
    return p.beta * (std::sqrt(A) - std::sqrt(p.A0)) / p.A0;
}

Real LinearElasticTubeLaw::pressureDerivative(Real A, const VesselParameters& p) const {
    return p.beta / (2.0 * p.A0 * std::sqrt(A));
}

Real LinearElasticTubeLaw::pressureFluxIntegral(Real A, const VesselParameters& p, Real rho) const {
    return p.beta / (3.0 * rho * p.A0) * (std::pow(A, 1.5) - std::pow(p.A0, 1.5));
}

Real LinearElasticTubeLaw::areaFromPressure(Real P, const VesselParameters& p) const {
    const Real sqrtA = (P * p.A0 / p.beta) + std::sqrt(p.A0);
    return sqrtA * sqrtA;
}

} // namespace hemo1d::physics