#pragma once

#include "hemo1d/core/network.hpp"
#include "hemo1d/core/vessel.hpp"
#include "hemo1d/physics/characteristics.hpp"
#include "hemo1d/physics/compatibility.hpp"
#include "hemo1d/physics/section_state.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {


// The 1D blood flow conservation law together with the wall constitutive model and
// the fluid properties. 
class BloodFlowModel {
public:
    BloodFlowModel(const TubeLaw& tubeLaw, FluidProperties fluid) noexcept
        : tubeLaw_(&tubeLaw), fluid_(fluid) {}
    
    // Wall and constitutive model functions

    // Transmural pressure P(A) - P0 for a vessel with the given parameters
    Real pressure(Real A, const VesselParameters& p) const;

    // Wave speed c(A) = sqrt((A / rho) dP/dA)
    Real waveSpeed(Real A, const VesselParameters& p) const;

    // Functions from hyperbolic system

    // Physical flux from the hyperbolic system F(U) = (Q, alpha Q^2 / A + pressureFluxIntegral(A))^T.
    SectionState physicalFlux(SectionState u, const VesselParameters& p) const;

    // Eigenvalues lambda = alpha u +/- c_alpha from the flux Jacobian, evaluated at U.
    Characteristics characteristics(SectionState u, const VesselParameters& p) const;

    // Left eigenvectors of the flux Jacobian, in the lambdaMinus, lambdaPlus order.
    LeftEigenvectors leftEigenvectors(SectionState u, const VesselParameters& p) const;

    // Forward Euler prediction of U at the next step from the local quasi-linear PDE:
    // dU/dt = -H(U) dU/dz - B(U).
    SectionState compatibilityPrediction(
        SectionState u, SectionGradient g, const VesselParameters& p, Real dt
    ) const;

    // Total (static + dynamic) pressure P(A) + 0.5 rho u^2.
    Real totalPressure(SectionState u, const VesselParameters& p) const;

    // Getters

    const TubeLaw& tubeLaw() const noexcept { return *tubeLaw_; }
    const FluidProperties& fluid() const noexcept { return fluid_; }
    Real density() const noexcept { return fluid_.density; }

private:
    const TubeLaw* tubeLaw_;
    FluidProperties fluid_;
};

} // namespace hemo1d::physics
