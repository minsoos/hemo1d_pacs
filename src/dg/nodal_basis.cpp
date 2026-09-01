#include "hemo1d/dg/nodal_basis.hpp"

#include <stdexcept>
#include <utility>

namespace hemo1d::dg {

NodalBasis::NodalBasis(std::vector<Real> nodes) : nodes_(std::move(nodes)) {
    if (nodes_.empty()) {
        throw std::invalid_argument("NodalBasis: at least one node is required");
    }
    computeBarycentricWeights();
    computeDifferentiationMatrix();
}

void NodalBasis::computeBarycentricWeights() {
    const Index n = nodes_.size();
    barycentricWeights_.assign(n, 1.0);
    for (Index j = 0; j < n; ++j) {
        for (Index k = 0; k < n; ++k) {
            if (k == j) continue;
            const Real diff = nodes_[j] - nodes_[k];
            if (diff == 0.0) {
                throw std::invalid_argument("NodalBasis: duplicate nodes are not allowed");
            }
            barycentricWeights_[j] /= diff;
        }
    }
}

void NodalBasis::computeDifferentiationMatrix() {
    const Index n = nodes_.size();
    diffMatrix_ = DenseMatrix::Zero(n, n);
    for (Index i = 0; i < n; ++i) {
        Real rowSum = 0.0;
        for (Index j = 0; j < n; ++j) {
            if (i == j) continue;
            const Real value = (
                (barycentricWeights_[j] / barycentricWeights_[i]) /
                (nodes_[i] - nodes_[j])
            );
            diffMatrix_(i, j) = value;
            rowSum += value;
        }
        diffMatrix_(i, i) = -rowSum;
    }
}


std::vector<Real> NodalBasis::evaluate(Real x) const {
    const Index n = nodes_.size();
    std::vector<Real> l(n, 0.0);

    for (Index i = 0; i < n; ++i) {
        if (x == nodes_[i]) {
            l[i] = 1.0;
            return l;
        }
    }

    std::vector<Real> terms(n);
    Real denom = 0.0;
    for (Index j = 0; j < n; ++j) {
        terms[j] = barycentricWeights_[j] / (x - nodes_[j]);
        denom += terms[j];
    }
    for (Index i = 0; i < n; ++i) {
        l[i] = terms[i] / denom;
    }

    return l;
}


std::vector<Real> NodalBasis::evaluateDerivative(Real x) const {
    const Index n = nodes_.size();

    for (Index i = 0; i < n; ++i) {
        if (x == nodes_[i]) {
            std::vector<Real> d(n);
            for (Index j = 0; j < n; ++j) d[j] = diffMatrix_(i, j);
            return d;
        }
    }

    const std::vector<Real> l = evaluate(x);
    Real s = 0.0;
    for (Index j = 0; j < n; ++j) {
        s += l[j] / (x - nodes_[j]);
    }

    std::vector<Real> d(n);
    for (Index i = 0; i < n; ++i) {
        d[i] = l[i] * (s - 1.0 / (x - nodes_[i]));
    }

    return d;
}

} // namespace hemo1d::dg