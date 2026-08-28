#include "hemo1d/physics/slope_limiter.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <span>
#include <vector>

#include "hemo1d/core/parallel.hpp"

namespace hemo1d::physics {

namespace {

Real minmod(std::span<const Real> values) {
    if (values.empty()) return 0.0;

    const Real sign = values.front() > 0.0 ? 1.0 : (values.front() < 0.0 ? -1.0 : 0.0);
    if (sign == 0.0) return 0.0;

    Real minAbs = std::abs(values.front());
    for (Real v : values) {
        if ((v > 0.0) != (sign > 0.0)) return 0.0;
        minAbs = std::min(minAbs, std::abs(v));
    }
    return sign * minAbs;
}

Real cellAverage(const dg::Element& el, const std::vector<Real>& component) {
    const dg::ReferenceElement& ref = el.referenceElement();
    const dg::QuadratureRule& quad = ref.quadrature();
    const DenseMatrix& L = ref.basisAtQuadrature();
    const Index n = el.numDofs();
    const Index offset = el.dofOffset();

    Real sum = 0.0;
    for (Index q = 0; q < quad.points.size(); ++q) {
        Real value = 0.0;
        for (Index i = 0; i < n; ++i) value += L(q, i) * component[offset + i];
        sum += quad.weights[q] * value;
    }
    return 0.5 * sum;
}

void limitComponent(
    const dg::Element& el, Real ownAvg, std::optional<Real> leftAvg,
    std::optional<Real> rightAvg, Real tolerance, std::vector<Real>& component
) {
    const Index n = el.numDofs();
    const Index offset = el.dofOffset();
    const Real uLeft = component[offset];
    const Real uRight = component[offset + n - 1];
    const Real h = el.length();

    std::array<Real, 3> candidates;
    Index count = 0;
    candidates[count++] = (uRight - uLeft) / h;
    if(rightAvg.has_value()) candidates[count++] = (*rightAvg - ownAvg) / h;
    if(leftAvg.has_value()) candidates[count++] = (ownAvg - *leftAvg) / h;

    const Real limitedSlope = minmod(std::span<const Real>(candidates.data(), count));
    const Real vLeft = ownAvg - limitedSlope * 0.5 * h;
    const Real vRight = ownAvg + limitedSlope * 0.5 * h;

    if (std::abs(vLeft - uLeft) < tolerance && std::abs(vRight - uRight) < tolerance) {
        return;
    }

    const Real xCenter = 0.5 * (el.zLeft() + el.zRight());
    const std::vector<Real>& nodes = el.physicalNodes();
    for (Index i = 0; i < n; ++i) {
        component[offset + i] = ownAvg + limitedSlope * (nodes[i] - xCenter);
    }
}

} // namespace


void MinmodLimiter::apply(State& state, const dg::Mesh& mesh) const {
    const std::vector<dg::Element>& elements = mesh.elements();
    const Index numElements = elements.size();

    std::vector<Real> avgA(numElements);
    std::vector<Real> avgQ(numElements);

# pragma omp parallel for if (numElements >= kOmpParallelThresholdLight)
    for (Index e = 0; e < numElements; ++e) {
        avgA[e] = cellAverage(elements[e], state.A);
        avgQ[e] = cellAverage(elements[e], state.Q);
    }

#pragma omp parallel for if (numElements >= kOmpParallelThresholdLight)
    for (Index e = 0; e < numElements; ++e) {
        const dg::Element& el = elements[e];
        const std::optional<Index> leftIdx = el.leftNeighbor();
        const std::optional<Index> rightIdx = el.rightNeighbor();

        const std::optional<Real> leftAvgA = leftIdx ? std::optional<Real>(avgA[*leftIdx]) : std::nullopt;
        const std::optional<Real> rightAvgA = rightIdx ? std::optional<Real>(avgA[*rightIdx]) : std::nullopt;
        limitComponent(el, avgA[e], leftAvgA, rightAvgA, tolerance_, state.A);

        const std::optional<Real> leftAvgQ = leftIdx ? std::optional<Real>(avgQ[*leftIdx]) : std::nullopt;
        const std::optional<Real> rightAvgQ = rightIdx ? std::optional<Real>(avgQ[*rightIdx]) : std::nullopt;
        limitComponent(el, avgQ[e], leftAvgQ, rightAvgQ, tolerance_, state.Q);
    }
}


} // namespace hemo1d::physics
