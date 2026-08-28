#pragma once

#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/state.hpp"

namespace hemo1d::physics {


// Stabilizes the DG solution by damping spurious oscillations
// near steep gradients.    
class SlopeLimiter {
public:
    virtual ~SlopeLimiter() = default;
    virtual void apply(State& state, const dg::Mesh& mesh) const = 0;
};


// Implements the minmod slope limiter, which is a simple and widely used limiter in DG methods.
class MinmodLimiter : public SlopeLimiter {
public:
    explicit MinmodLimiter(Real tolerance = 1e-8) : tolerance_(tolerance) {};

    void apply(State& state, const dg::Mesh& mesh) const override;

private:
    Real tolerance_;
};


} // namespace hemo1d::physics
