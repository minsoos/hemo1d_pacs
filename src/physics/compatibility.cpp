#include "hemo1d/physics/compatibility.hpp"
 
#include "hemo1d/physics/characteristics.hpp"
 
namespace hemo1d::physics {
 
LeftEigenvectors computeLeftEigenvectorscomputeLeftEigenvectors(Real A, Real Q, Real A0,
                    Real beta, Real alpha, Real rho, const TubeLaw& tubelaw){
    const Real u = Q / A;
    const Characteristics c = computeCharacteristics(A, Q, A0, beta, alpha, rho, tubelaw);
    const Real cAlpha = 0.5 * (c.lambdaPlus - c.lambdaMinus);
    return {{cAlpha - alpha * u, 1.0}, {-cAlpha - alpha * u, 1.0}};
}

std::pair<Real, Real> compatibilityPrediction(Real A, Real Q, Real dAdz, Real dQdz, 
                    Real A0, Real beta, Real alpha, Real rho, Real frictionKr, 
                    const TubeLaw& tubelaw, Real dt){

                    const Real u = Q / A;
                    const Real c = tubelaw.waveSpeed(A, A0, beta, rho);
                    const Real hA = dQdz;
                    const Real hQ = (c * c - alpha * u * u) * dAdz + 2.0 * alpha * u * dQdz;
                    const Real bQ = frictionKr * u;
                    return {A - dt * hA, Q - dt * (hQ + bQ)};   


}

} // namespace hemo1d::physics