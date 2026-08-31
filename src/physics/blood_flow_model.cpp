#include "hemo1d/physics/blood_flow_model.hpp"

#include <cmath>

#include "hemo1d/physics/conservation_law.hpp"

namespace hemo1d::physics {

Real BloodFlowModel::pressure(Real A, const VesselParameters& p) const {
    return tubeLaw_->pressure(A, p.A0, p.beta);
}

Real BloodFlowModel::waveSpeed(Real A, const VesselParameters& p) const {
    return tubeLaw_->waveSpeed(A, p.A0, p.beta, fluid_.density);
}

SectionState BloodFlowModel::physicalFlux(SectionState u, const VesselParameters& p) const {
    return {u.Q, p.alpha * u.Q * u.Q / u.A + tubeLaw_->pressureFluxIntegral(u.A, p.A0, p.beta, fluid_.density)};
}

Characteristics BloodFlowModel::characteristics(SectionState u, const VesselParameters& p) const {
    const Real vel = u.velocity();
    const Real c = tubeLaw_->waveSpeed(u.A, p.A0, p.beta, fluid_.density);
    const Real cAlpha = std::sqrt(c * c + p.alpha * (p.alpha - 1.0) * vel * vel);

    return {p.alpha * vel - cAlpha, p.alpha * vel + cAlpha};
}

LeftEigenvectors BloodFlowModel::leftEigenvectors(SectionState u, const VesselParameters& p) const {
    const Real vel = u.velocity();
    const Characteristics c = characteristics(u, p);
    const Real cAlpha = 0.5 * (c.lambdaPlus - c.lambdaMinus);
    return {{cAlpha - p.alpha * vel, 1.0}, {-cAlpha - p.alpha * vel, 1.0}};
}

SectionState BloodFlowModel::compatibilityPrediction(
    SectionState u, SectionGradient g, const VesselParameters& p, Real dt
) const {
    const Real vel = u.velocity();
    const Real c = tubeLaw_->waveSpeed(u.A, p.A0, p.beta, fluid_.density);
    const Real hA = g.dQdz;
    const Real hQ = (c * c - p.alpha * vel * vel) * g.dAdz + 2.0 * p.alpha * vel * g.dQdz;
    const Real bQ = p.frictionKr * vel;
    return {u.A - dt * hA, u.Q - dt * (hQ + bQ)};
}

Real BloodFlowModel::totalPressure(SectionState u, const VesselParameters& p) const {
    const Real vel = u.velocity();
    return tubeLaw_->pressure(u.A, p.A0, p.beta) + 0.5 * fluid_.density * vel * vel;
}

} // namespace hemo1d::physics
