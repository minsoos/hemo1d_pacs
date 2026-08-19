#pragma once

#include <vector>

#include "hemo1d/core/types.hpp"

namespace hemo1d {

class DenseMatrix {
public:
    DenseMatrix() = default;
    DenseMatrix(Index rows, Index cols, Real fill = 0.0)
        : rows_(rows), cols_(cols), data_(rows * cols, fill) {}

    Index rows() const noexcept { return rows_; }
    Index cols() const noexcept { return cols_; }

    Real& operator()(Index row, Index col) { return data_[row * cols_ + col]; }
    Real operator()(Index row, Index col) const { return data_[row * cols_ + col]; }

private:
    Index rows_ = 0;
    Index cols_ = 0;
    std::vector<Real> data_;
};

} // namespace hemo1d