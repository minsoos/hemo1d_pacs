#pragma once

#include <vector>

#include "hemo1d/core/dense_matrix.hpp"

namespace hemo1d {

std::vector<Real> solveLinearSystem(
    DenseMatrix A, std::vector<Real> b
);

DenseMatrix invertDense(const DenseMatrix& A);

} // namespace hemo1d