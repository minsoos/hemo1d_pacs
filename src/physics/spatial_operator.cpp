#include "hemo1d/physics/spatial_operator.hpp"

#include <array>
#include <stdexcept>
#include <vector>

#include "hemo1d/core/parallel.hpp"
#include "hemo1d/physics/conservation_law.hpp"

namespace hemo1d::physics {

SpatialOperator::SpatialOperator(
    const dg::Mesh& mesh, FluidProperties fluid,
    const TubeLaw& tubeLaw, const NumericalFlux& flux,
    const BoundaryStateProvider& boundaryProvider
)
    : mesh_(mesh),
    fluid_(fluid),
    tubeLaw_(tubeLaw),
    flux_(flux),
    boundaryProvider_(boundaryProvider) 
{
    for (const dg::Element& el : mesh_.elements()) {
        if (el.numDofs() > kMaxElementDofs) {
            throw std::runtime_error(
                "SpatialOperator: element polynomial order exceeds the supported maximum (hemo1d::kMaxElementDofs)"
            );
        }
    }
}

void SpatialOperator::evaluate(const State& u, Real time, State& dudt) const {
    const std::vector<dg::Element>& elements = mesh_.elements();
    const Index numElements = elements.size();

#pragma omp parallel for if (numElements >= kOmpParallelThreshold)
    for (Index e = 0; e < numElements; ++e) {
        const dg::Element& el = elements[e];
        const VesselParameters& p = el.vesselParameters();
        const Index n = el.numDofs();
        const Index offset = el.dofOffset();
        const Real J = el.jacobian();

        const Real Al = u.A[offset];
        const Real Ql = u.Q[offset];
        const Real Ar = u.A[offset + n - 1];
        const Real Qr = u.Q[offset + n - 1];

        std::pair<Real, Real> outsideLeft;
        if (el.leftNeighbor().has_value()) {
            const dg::Element& neighbor = elements[*el.leftNeighbor()];
            const Index nOff = neighbor.dofOffset() + neighbor.numDofs() - 1;
            outsideLeft = {u.A[nOff], u.Q[nOff]};
        } else {
            outsideLeft = boundaryProvider_.exteriorState(el.vesselId(), VesselEnd::Proximal, time);
        }

        std::pair<Real, Real> outsideRight;
        if (el.rightNeighbor().has_value()) {
            const dg::Element& neighbor = elements[*el.rightNeighbor()];
            const Index nOff = neighbor.dofOffset();
            outsideRight = { u.A[nOff], u.Q[nOff] };
        } else {
            outsideRight = boundaryProvider_.exteriorState(el.vesselId(), VesselEnd::Distal, time);
        }

        const auto [leftFluxA, leftFluxQ] = flux_.compute(
            outsideLeft.first, outsideLeft.second, Al, Ql, p.A0, p.beta, p.alpha, fluid_.density, tubeLaw_);
        const auto [rightFluxA, rightFluxQ] = flux_.compute(
            Ar, Qr, outsideRight.first, outsideRight.second, p.A0, p.beta, p.alpha, fluid_.density, tubeLaw_);

        const dg::ReferenceElement& ref = el.referenceElement();
        const dg::QuadratureRule& quad = ref.quadrature();
        const DenseMatrix& L = ref.basisAtQuadrature();
        const DenseMatrix& dL = ref.basisDerivativeAtQuadrature();
        const Index nq = quad.points.size();

        std::array<Real, kMaxElementDofs> rhsA{};
        std::array<Real, kMaxElementDofs> rhsQ{};

        for (Index q = 0; q < nq; ++q) {
            Real Aq = 0.0;
            Real Qq = 0.0;
            for (Index i = 0; i < n; ++i) {
                Aq += L(q, i) * u.A[offset + i];
                Qq += L(q, i) * u.Q[offset + i];
            }
            const auto [fA, fQ] = physicalFlux(Aq, Qq, p.A0, p.beta, p.alpha, fluid_.density, tubeLaw_);
            const Real sourceQ = -p.frictionKr * Qq / Aq;
            const Real w = quad.weights[q];
            for (Index i = 0; i < n; ++i) {
                rhsA[i] += w * fA * dL(q, i);
                rhsQ[i] += w * fQ * dL(q, i) + J * w * sourceQ * L(q, i);
            }
        }

        rhsA[0] += leftFluxA;
        rhsA[n - 1] -= rightFluxA;
        rhsQ[0] += leftFluxQ;
        rhsQ[n - 1] -= rightFluxQ;

        const DenseMatrix& Minv = ref.massMatrixInverse();
        for (Index i = 0; i < n; ++i) {
            Real dA = 0.0;
            Real dQ = 0.0;
            for (Index j = 0; j < n; ++j) {
                dA += Minv(i, j) * rhsA[j];
                dQ += Minv(i, j) * rhsQ[j];
            }
            dudt.A[offset + i] = dA / J;
            dudt.Q[offset + i] = dQ / J;
        }
    }
}

} // namespace hemo1d::physics
