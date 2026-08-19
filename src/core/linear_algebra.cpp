#include "hemo1d/core/linear_algebra.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace hemo1d {

std::vector<Real> solveLinearSystem(
    DenseMatrix A, std::vector<Real> b
) {
    const Index n = A.rows();
    if (A.cols() != n || b.size() != n) {
        throw std::invalid_argument("solveLinearSystem: dimension mismatch");
    }
    
    for (Index col = 0; col < n; ++col) {
        const Real pivot = A(col, col);
        for (Index row = col + 1; row < n; ++row) {
            const Real factor = A(row, col) / pivot;
            for (Index k = col; k < n; ++k) A(row, k) -= factor * A(col, k);
            b[row] -= factor * b[col];
        }
    }

    std::vector<Real> x(n);
    for (Index i = n; i-- > 0;) {
        Real sum = b[i];
        for (Index k = i + 1; k < n; ++k) sum -= A(i, k) * x[k];
        x[i] = sum / A(i, i);
    }

    return x;
}


DenseMatrix invertDense(const DenseMatrix& A) {
    const Index n = A.rows();
    DenseMatrix inv(n, n);
    for (Index col = 0; col < n; ++col) {
        std::vector<Real> unit(n, 0.0);
        unit[col] = 1.0;
        const std::vector<Real> x = solveLinearSystem(A, unit);
        for (Index row = 0; row < n; ++row) inv(row, col) = x[row];
    }
    return inv;
}

} // namespace hemo1d
