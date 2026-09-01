#pragma once

#include <Eigen/Dense>

#include "hemo1d/core/types.hpp"

namespace hemo1d {

using DenseMatrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

} // namespace hemo1d