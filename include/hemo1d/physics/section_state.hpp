#pragma once

#include "hemo1d/core/types.hpp"

namespace hemo1d::physics {

struct SectionState {
    Real A = 0.0;
    Real Q = 0.0;

    constexpr Real velocity() const { return Q / A; }
};

struct SectionGradient {
    Real dAdz = 0.0;
    Real dQdz = 0.0;
};

} // namespace hemo1d::physics
