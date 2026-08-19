#include "hemo1d/physics/flux.hpp"

#include <algorithm>
#include <cmath>

#include "hemo1d/physics/conservation_law.hpp"
#include "hemo1d/physics/characteristics.hpp"

namespace hemo1d::physics {

std::pair<Real, Real> LaxFriedrichsFlux::compute(
    Real AL, Real QL, Real AR, Real QR,
    Real A0, Real beta, Real alpha, Real rho, const TubeLaw& tubeLaw
) const {

    const auto [FAL, FQL] = physicalFlux(AL, QL, A0, beta, alpha, rho, tubeLaw);
    const auto [FAR, FQR] = physicalFlux(AR, QR, A0, beta, alpha, rho, tubeLaw);

    const Characteristics leftChar = computeCharacteristics(AL, QL, A0, beta, alpha, rho, tubeLaw);
    const Characteristics rightChar = computeCharacteristics(AR, QR, A0, beta, alpha, rho, tubeLaw);
    const Real lambdaMax = std::max({std::abs(leftChar.lambdaMinus), std::abs(leftChar.lambdaPlus),
                                    std::abs(rightChar.lambdaMinus), std::abs(rightChar.lambdaPlus)});

    return {
        0.5 * (FAL + FAR) - 0.5 * lambdaMax * (AR - AL),
        0.5 * (FQL + FQR) - 0.5 * lambdaMax * (QR - QL)
    };
}


std::pair<Real, Real> HllFlux::compute(
    Real AL, Real QL, Real AR, Real QR,
    Real A0, Real beta, Real alpha, Real rho, const TubeLaw& tube
) const {

    const Characteristics leftChar = computeCharacteristics(AL, QL, A0, beta, alpha, rho, tube);
    const Characteristics rightChar = computeCharacteristics(AR, QR, A0, beta, alpha, rho, tube);

    const Real sL = std::min(leftChar.lambdaMinus, rightChar.lambdaMinus);
    const Real sR = std::max(leftChar.lambdaPlus, rightChar.lambdaPlus);

    const auto [FAL, FQL] = physicalFlux(AL, QL, A0, beta, alpha, rho, tube);
    if (sL >= 0.0) return {FAL, FQL};

    const auto [FAR, FQR] = physicalFlux(AR, QR, A0, beta, alpha, rho, tube);
    if (sR <= 0.0) return {FAR, FQR};

    const Real denom = sR - sL;

    return {
        (sR * FAL - sL * FAR + sR * sL * (AR - AL)) / denom,
        (sR * FQL - sL * FQR + sR * sL * (QR - QL)) / denom
    };

}

} // namespace hemo1d::physics