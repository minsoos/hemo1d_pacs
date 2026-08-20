#pragma once

#include <vector>

#include "hemo1d/core/types.hpp"

namespace hemo1d::dg {

struct QuadratureRule {
    std::vector<Real> points;
    std::vector<Real> weights;
};

QuadratureRule gaussLegendre(unsigned numPoints);

QuadratureRule gaussLobattoLegendre(unsigned numPoints);

} // namespace hemo1d::dg