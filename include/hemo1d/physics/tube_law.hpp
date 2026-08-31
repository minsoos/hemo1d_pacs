#pragma once

#include "hemo1d/core/types.hpp"
#include "hemo1d/core/vessel.hpp"

namespace hemo1d::physics {

// Abstracts the constitutive relation between transmural pressure and
// cross-sectional area (the "tube law"), so alternative wall models can be
// substituted.
class TubeLaw {
public:
    virtual ~TubeLaw() = default;

    virtual Real pressure(Real A, const VesselParameters& p) const = 0;
    virtual Real pressureDerivative(Real A, const VesselParameters& p) const = 0;

    virtual Real pressureFluxIntegral(Real A, const VesselParameters& p, Real rho) const = 0;
    
    virtual Real areaFromPressure(Real P, const VesselParameters& p) const = 0;

    Real waveSpeed(Real A, const VesselParameters& p, Real rho) const;
};


// The linear elastic (membrane) tube law:
// P - P0 = beta * (sqrt(A) - sqrt(A0)) / A0
class LinearElasticTubeLaw : public TubeLaw {
public:
    Real pressure(Real A, const VesselParameters& p) const override;
    Real pressureDerivative(Real A, const VesselParameters& p) const override;
    Real pressureFluxIntegral(Real A, const VesselParameters& p, Real rho) const override;
    Real areaFromPressure(Real P, const VesselParameters& p) const override;
};

}  // namespace hemo1d::physics