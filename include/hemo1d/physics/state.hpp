#pragma once
 
#include <vector>
 
#include "hemo1d/core/types.hpp"
 
namespace hemo1d::physics {
 
// State::A[i], State::Q[i] for obtaining the i-th element in the state
struct State {
    std::vector<Real> A;
    std::vector<Real> Q;
 
    State() = default;
    explicit State(Index totalDofs) : A(totalDofs, 0.0), Q(totalDofs, 0.0) {}
 
    Index size() const noexcept { return A.size(); }
};
 
} // namespace hemo1d::physics
 