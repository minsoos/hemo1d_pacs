#pragma once

#include <utility>

#include "hemo1d/core/types.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {

class NumericalFlux {
public:
    virtual ~NumericalFlux() = default;

    virtual std::pair<Real, Real> compute(
        Real AL, Real QL, Real AR, Real QR, 
        Real A0, Real beta, Real alpha, Real rho, const TubeLaw& tubeLaw
    ) const = 0;
};


class LaxFriedrichsFlux : public NumericalFlux {
public:
    std::pair<Real, Real> compute(
        Real AL, Real QL, Real AR, Real QR,
        Real A0, Real beta, Real alpha, Real rho, const TubeLaw& tubeLaw
    ) const override; 
};


class HllFlux : public NumericalFlux {
public:
    std::pair<Real, Real> compute(
        Real AL, Real QL, Real AR, Real QR,
        Real A0, Real beta, Real alpha, Real rho, const TubeLaw& tube
    ) const override;
};

} // namespace hemo1d::physics