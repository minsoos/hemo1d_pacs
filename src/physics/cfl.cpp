#include "hemo1d/physics/cfl.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "hemo1d/core/parallel.hpp"
#include "hemo1d/physics/characteristics.hpp"

namespace hemo1d::physics {

Real maxWaveSpeed(
    const State& u, const dg::Mesh& mesh, const FluidProperties& fluid,
    const TubeLaw& tubeLaw
) {
    const std::vector<dg::Element>& elements = mesh.elements();
    const Index numElements = elements.size();
    Real result = 0.0;

    #pragma omp parallel for reduction(max: result)
    for (Index e = 0; e < numElements; ++e) {
        const dg::Element& el = elements[e];
        const VesselParameters& p = el.vesselParameters();

        for (Index i = 0; i < el.numDofs(); ++i) {
            const Index idx = el.dofOffset() + i;
            const Characteristics c = computeCharacteristics(
                u.A[idx], u.Q[idx], p.A0, p.beta, p.alpha, fluid.density, tubeLaw
            );

            result = std::max({result, std::abs(c.lambdaMinus), std::abs(c.lambdaPlus)});
        }
    }
    return result;
}


Real cflTimeStep(
    const State& u, const dg::Mesh& mesh, const FluidProperties& fluid,
    const TubeLaw& tubeLaw, Real cflNumber
) {
    const std::vector<dg::Element>& elements = mesh.elements();
    const Index numElements = elements.size();
    Real minDt = std::numeric_limits<Real>::max();

    #pragma omp parallel for reduction(min: minDt)
    for (Index e = 0; e < numElements; ++e) {
        const dg::Element& el = elements[e];
        const VesselParameters& p = el.vesselParameters();
        Real localMaxSpeed = 0.0;

        for (Index i = 0; i < el.numDofs(); ++i) {
            const Index idx = el.dofOffset() + i;
            const Characteristics c = computeCharacteristics(
                u.A[idx], u.Q[idx], p.A0, p.beta, p.alpha, fluid.density, tubeLaw
            );
            localMaxSpeed = std::max({localMaxSpeed, std::abs(c.lambdaMinus), std::abs(c.lambdaPlus)});
        }
        if (localMaxSpeed > 0.0) {
            const Real bound = el.length() / ((2.0 * el.order() + 1.0) * localMaxSpeed);
            minDt = std::min(minDt, bound);
        }
    }
    return cflNumber * minDt;
}

} // namespace hemo1d::physics