#pragma once

#include <vector>

#include "hemo1d/core/dense_matrix.hpp"
#include "hemo1d/core/types.hpp"

namespace hemo1d::dg {

class NodalBasis {
public:
    explicit NodalBasis(std::vector<Real> nodes);

    unsigned order() const noexcept { return static_cast<unsigned>(nodes_.size() - 1); }
    Index numNodes() const noexcept { return nodes_.size(); }
    const std::vector<Real>& nodes() const noexcept { return nodes_; }

    std::vector<Real> evaluate(Real x) const;

    std::vector<Real> evaluateDerivative(Real x) const;

    const DenseMatrix& differentiationMatrix() const noexcept { return diffMatrix_; }

private:
    void computeBarycentricWeights();
    void computeDifferentiationMatrix();

    std::vector<Real> nodes_;
    std::vector<Real> barycentricWeights_;
    DenseMatrix diffMatrix_;
};

} // namespace hemo1d::dg