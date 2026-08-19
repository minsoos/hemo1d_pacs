#pragma once

#include "hemo1d/core/types.hpp"

namespace hemo1d::physics {

// Abstracts the constitutive relation between transmural pressure and
// cross-sectional area (the "tube law"), so alternative wall models can be
// substituted.
class TubeLaw {
public:
    virtual ~TubeLaw() = default;

    virtual Real pressure(Real A, Real A0, Real beta) const = 0;
    virtual Real pressureDerivative(Real A, Real A0, Real beta) const = 0;

    virtual Real pressureFluxIntegral(Real A, Real A0, Real beta, Real rho) const = 0;
    
    virtual Real areaFromPressure(Real P, Real A0, Real beta) const = 0;

    Real waveSpeed(Real A, Real A0, Real beta, Real rho) const;
};


// The linear elastic (membrane) tube law:
// P - P0 = beta * (sqrt(A) - sqrt(A0)) / A0
class LinearElasticTubeLaw : public TubeLaw {
public:
    Real pressure(Real A, Real A0, Real beta) const override;
    Real pressureDerivative(Real A, Real A0, Real beta) const override;
    Real pressureFluxIntegral(Real A, Real A0, Real beta, Real rho) const override;
    Real areaFromPressure(Real P, Real A0, Real beta) const override;
};

}  // namespace hemo1d::physics