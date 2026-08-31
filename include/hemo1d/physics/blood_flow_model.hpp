#pragma once

#include <utility>
#include <cmath>

#include "hemo1d/core/network.hpp"
#include "hemo1d/core/vessel.hpp"
#include "hemo1d/physics/section_state.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {

// The two characteristic (eigenvalue) speeds of the 1D
// blood flow hyperbolic equations at a given state.
struct Characteristics {
    Real lambdaMinus;  // alpha * U - c_alpha
    Real lambdaPlus;   // alpha * U + c_alpha
};

struct LeftEigenvectors {
    std::pair<Real, Real> lPlus;
    std::pair<Real, Real> lMinus;
};

// The 1D blood flow conservation law together with the wall constitutive model and
// the fluid properties. 
class BloodFlowModel {
public:
    BloodFlowModel(const TubeLaw& tubeLaw, FluidProperties fluid) noexcept
        : tubeLaw_(&tubeLaw), fluid_(fluid) {}
    
    // Wall and constitutive model functions

    // Transmural pressure P(A) - P0 for a vessel with the given parameters
    Real pressure(Real A, const VesselParameters& p) const { return tubeLaw_->pressure(A, p); }

    // Wave speed c(A) = sqrt((A / rho) dP/dA)
    Real waveSpeed(Real A, const VesselParameters& p) const {
        return tubeLaw_->waveSpeed(A, p, fluid_.density);
    }

    // Functions from hyperbolic system

    // Physical flux from the hyperbolic system F(U) = (Q, alpha Q^2 / A + pressureFluxIntegral(A))^T.
    SectionState physicalFlux(SectionState u, const VesselParameters& p) const {
        return {
            u.Q,
            p.alpha * u.Q * u.Q / u.A + tubeLaw_->pressureFluxIntegral(u.A, p, fluid_.density)
        };
    }

    // Eigenvalues lambda = alpha u +/- c_alpha from the flux Jacobian, evaluated at U.
    Characteristics characteristics(SectionState u, const VesselParameters& p) const {
        const Real vel = u.velocity();
        const Real c = tubeLaw_->waveSpeed(u.A, p, fluid_.density);
        const Real cAlpha = std::sqrt(c * c + p.alpha * (p.alpha - 1.0) * vel * vel);
        return {p.alpha * vel - cAlpha, p.alpha * vel +cAlpha};
    }

    // Left eigenvectors of the flux Jacobian, in the lambdaMinus, lambdaPlus order.
    LeftEigenvectors leftEigenvectors(SectionState u, const VesselParameters& p) const {
        const Real vel = u.velocity();
        const Characteristics c = characteristics(u, p);
        const Real cAlpha = 0.5 * (c.lambdaPlus - c.lambdaMinus);
        return {{cAlpha - p.alpha * vel, 1.0}, {-cAlpha - p.alpha * vel, 1.0}};
    }

    // Forward Euler prediction of U at the next step from the local quasi-linear PDE:
    // dU/dt = -H(U) dU/dz - B(U).
    SectionState compatibilityPrediction(
        SectionState u, SectionGradient g, const VesselParameters& p, Real dt
    ) const {
        const Real vel = u.velocity();
        const Real c = tubeLaw_->waveSpeed(u.A, p, fluid_.density);
        const Real hA = g.dQdz;
        const Real hQ = (c * c - p.alpha * vel * vel) * g.dAdz + 2.0 * p.alpha * vel * g.dQdz;
        const Real bQ = p.frictionKr * vel;
        return {u.A - dt * hA, u.Q - dt * (hQ + bQ)};
    }

    // Total (static + dynamic) pressure P(A) + 0.5 rho u^2.
    Real totalPressure(SectionState u, const VesselParameters& p) const {
        const Real vel = u.velocity();
        return tubeLaw_->pressure(u.A, p) + 0.5 * fluid_.density * vel * vel;
    }

    // Getters

    const TubeLaw& tubeLaw() const noexcept { return *tubeLaw_; }
    const FluidProperties& fluid() const noexcept { return fluid_; }
    Real density() const noexcept { return fluid_.density; }

private:
    const TubeLaw* tubeLaw_;
    FluidProperties fluid_;
};

} // namespace hemo1d::physics
