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
        Index pivotRow = col;
        Real pivotMag = std::abs(A(col, col));
        for (Index row = col + 1; row < n; ++row) {
            if (std::abs(A(row, col)) > pivotMag) {
                pivotMag = std::abs(A(row, col));
                pivotRow = row;
            }
        }
        if (pivotMag < 1e-14) {
            throw std::runtime_error("solveLinearSystem: matrix is (numerically) singular");
        }
        if (pivotRow != col) {
            for (Index k = 0; k < n; ++k) std::swap(A(col, k), A(pivotRow, k));
            std::swap(b[col], b[pivotRow]);
        }
        for (Index row = col + 1; row < n; ++row) {
            const Real factor = A(row, col) / A(col, col);
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
