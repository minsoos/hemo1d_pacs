#include "hemo1d/dg/reference_element.hpp"


namespace hemo1d::dg {

namespace {

NodalBasis makeGllBasis(unsigned order) { 
    return NodalBasis(gaussLobattoLegendre(order + 1).points);
}

} // namespace

ReferenceElement::ReferenceElement(unsigned order)
    : basis_(makeGllBasis(order)), quadrature_(gaussLegendre(order + 1)) {
    
    const Index n = basis_.numNodes();
    const Index nq = quadrature_.points.size();

    basisAtQuadrature_ = DenseMatrix::Zero(nq, n);
    basisDerivativeAtQuadrature_ = DenseMatrix::Zero(nq, n);
    for (Index q = 0; q < nq; ++q) {
        const std::vector<Real> l = basis_.evaluate(quadrature_.points[q]);
        const std::vector<Real> dl = basis_.evaluateDerivative(quadrature_.points[q]);
        for (Index i = 0; i < n; ++i) {
            basisAtQuadrature_(q, i) = l[i];
            basisDerivativeAtQuadrature_(q, i) = dl[i];
        }
    }

    massMatrix_ = DenseMatrix::Zero(n, n);
    for (Index i = 0; i < n; ++i) {
        for (Index j = 0; j < n; ++j) {
            Real sum = 0.0;
            for (Index q = 0; q < nq; ++q) {
                sum += quadrature_.weights[q] * basisAtQuadrature_(q, i) * basisAtQuadrature_(q, j);
            }
            massMatrix_(i, j) = sum;
        }
    }

    massMatrixInverse_ = massMatrix_.partialPivLu().solve(DenseMatrix::Identity(n, n));
}


} // namespace hemo1d::dg