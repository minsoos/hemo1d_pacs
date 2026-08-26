#pragma once

#include "hemo1d/core/dense_matrix.hpp"
#include "hemo1d/dg/nodal_basis.hpp"
#include "hemo1d/dg/quadrature.hpp"

namespace hemo1d::dg {

class ReferenceElement {
public:
    explicit ReferenceElement(unsigned order);

    unsigned order() const noexcept { return basis_.order(); }
    Index numNodes() const noexcept { return basis_.numNodes(); }
    const std::vector<Real>& nodes() const noexcept { return basis_.nodes(); }
    const NodalBasis& basis() const noexcept { return basis_; }

    const DenseMatrix& differentiationMatrix() const noexcept { 
        return basis_.differentiationMatrix(); 
    }

    const DenseMatrix& massMatrix() const noexcept { return massMatrix_; }

    const DenseMatrix& massMatrixInverse() const noexcept { return massMatrixInverse_; }

    const QuadratureRule& quadrature() const noexcept { return quadrature_; }

    const DenseMatrix& basisAtQuadrature() const noexcept { return basisAtQuadrature_; }
    const DenseMatrix& basisDerivativeAtQuadrature() const noexcept { 
        return basisDerivativeAtQuadrature_; 
    }

private:
    NodalBasis basis_;
    QuadratureRule quadrature_;
    DenseMatrix massMatrix_;
    DenseMatrix massMatrixInverse_;
    DenseMatrix basisAtQuadrature_;
    DenseMatrix basisDerivativeAtQuadrature_;
};

} // namespace hemo1d::dg