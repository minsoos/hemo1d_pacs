#pragma once

#include "hemo1d/physics/boundary_state_provider.hpp"
#include "hemo1d/physics/flux.hpp"
#include "hemo1d/physics/slope_limiter.hpp"
#include "hemo1d/physics/spatial_operator.hpp"
#include "hemo1d/physics/state.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics {

class Solver {
public:
    Solver(
        const dg::Mesh& mesh, FluidProperties fluid, const TubeLaw& tubeLaw,
        const NumericalFlux& flux, const BoundaryStateProvider& boundaryProvider,
        const SlopeLimiter* limiter = nullptr
    );

    void step(State& u, Real time, Real dt) const;

private:
    SpatialOperator operatorL_;
    const SlopeLimiter* limiter_;
    mutable State stage_;
    mutable State k1_;
    mutable State k2_;
};

} // namespace hemo1d::physics
