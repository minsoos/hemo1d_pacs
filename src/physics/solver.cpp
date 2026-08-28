#include "hemo1d/physics/solver.hpp"

#include "hemo1d/core/parallel.hpp"

namespace hemo1d::physics {

Solver::Solver(
    const dg::Mesh& mesh, FluidProperties fluid, const TubeLaw& tubeLaw,
    const NumericalFlux& flux, const BoundaryStateProvider& boundaryProvider,
    const SlopeLimiter* limiter
) 
    : operatorL_(mesh, fluid, tubeLaw, flux, boundaryProvider),
    limiter_(limiter),
    stage_(mesh.totalDofs()),
    k1_(mesh.totalDofs()),
    k2_(mesh.totalDofs())
{}


void Solver::step(State& u, Real time, Real dt) const {
    const Index n = u.size();

    operatorL_.evaluate(u, time, k1_);

#pragma omp parallel for if (n >= kOmpParallelThreshold)
    for (Index i = 0; i < n; ++i) {
        stage_.A[i] = u.A[i] + dt * k1_.A[i];
        stage_.Q[i] = u.Q[i] + dt * k1_.Q[i];
    }
    if (limiter_) limiter_->apply(stage_, operatorL_.mesh());

    operatorL_.evaluate(stage_, time + dt, k2_);

#pragma omp parallel for if (n >= kOmpParallelThreshold)
    for (Index i = 0; i < n; ++i) {
        u.A[i] += 0.5 * dt * (k1_.A[i] + k2_.A[i]);
        u.Q[i] += 0.5 * dt * (k1_.Q[i] + k2_.Q[i]);
    }
    if (limiter_) limiter_->apply(u, operatorL_.mesh());
}

} // namespace hemo1d::physics
