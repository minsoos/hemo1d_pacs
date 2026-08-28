// Benchmarks Hemo1D's hot paths:
//   1. SpatialOperator::evaluate() alone, across a range of element counts
//      (used to calibrate hemo1d::kOmpParallelThreshold,
//      include/hemo1d/core/parallel.hpp).
//   2. A full Solver::step() on a realistic Y-bifurcation network, broken
//      down by component (SpatialOperator, MinmodLimiter, BoundarySolver,
//      the RK combine loops), to find which part of a real timestep
//      actually dominates wall-clock time.
//
// Usage: hemo1d_bench [repeats]
 
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
 
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/boundary_solver.hpp"
#include "hemo1d/physics/flux.hpp"
#include "hemo1d/physics/slope_limiter.hpp"
#include "hemo1d/physics/solver.hpp"
#include "hemo1d/physics/spatial_operator.hpp"
#include "hemo1d/physics/tube_law.hpp"
 
using namespace hemo1d;
using namespace hemo1d::dg;
using namespace hemo1d::physics;
 
namespace {
 
class ConstantBoundaryStateProvider : public BoundaryStateProvider {
public:
    ConstantBoundaryStateProvider(Real A, Real Q) : A_(A), Q_(Q) {}
    std::pair<Real, Real> exteriorState(Id, VesselEnd, Real) const override { return {A_, Q_}; }
 
private:
    Real A_;
    Real Q_;
};
 
double msSince(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}
 
// ---------------------------------------------------------------------
// Part 1: SpatialOperator::evaluate() alone, across element counts.
// ---------------------------------------------------------------------
 
Network makeSingleVesselNetwork(Index nElements) {
    VesselParameters params;
    params.length = static_cast<Real>(nElements); // keep element size ~1 regardless of count
    params.A0 = 0.126;
    params.beta = 6.06e5;
    params.alpha = 4.0 / 3.0;
    params.frictionKr = 0.0;
    params.nElements = nElements;
 
    std::vector<Vessel> vessels{Vessel(1, "v", params)};
    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, bc),
        Node(2, "out", {{1, VesselEnd::Distal}}, bc),
    };
    return Network(FluidProperties{}, std::move(vessels), std::move(nodes));
}
 
double benchmarkEvaluate(Index nElements, int repeats) {
    const Network network = makeSingleVesselNetwork(nElements);
    const Mesh mesh(network);
 
    State u(mesh.totalDofs());
    State dudt(mesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = 0.126;
        u.Q[i] = 0.05;
    }
 
    LinearElasticTubeLaw law;
    HllFlux flux;
    ConstantBoundaryStateProvider provider(0.126, 0.05);
    SpatialOperator op(mesh, FluidProperties{}, law, flux, provider);
 
    op.evaluate(u, 0.0, dudt); // warm-up (page faults, cache fill)
 
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) {
        op.evaluate(u, 0.0, dudt);
    }
    return msSince(start) / repeats;
}
 
void runEvaluateSweep(int repeats) {
    const std::vector<Index> sizes = {8,    16,   32,   64,    128,   256,  512,
                                       1024, 2048, 4096, 8192, 16384, 32768};
 
    std::printf("=== Part 1: SpatialOperator::evaluate() vs element count ===\n");
    std::printf("%10s %15s\n", "elements", "ms/evaluate");
    for (Index n : sizes) {
        const double ms = benchmarkEvaluate(n, repeats);
        std::printf("%10zu %15.6f\n", static_cast<size_t>(n), ms);
    }
    std::printf("\n");
}
 
// ---------------------------------------------------------------------
// Part 1b: MinmodLimiter::apply() alone, across element counts. Its
// per-element cost is much lower than SpatialOperator's (no quadrature
// loop, no small matrix solves), so it needs its own, separately tuned
// dispatch threshold (kOmpParallelThresholdLight) -- see parallel.hpp.
// ---------------------------------------------------------------------
 
double benchmarkLimiter(Index nElements, int repeats) {
    const Network network = makeSingleVesselNetwork(nElements);
    const Mesh mesh(network);
 
    State u(mesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = 0.126;
        u.Q[i] = 0.05;
    }
 
    MinmodLimiter limiter;
    limiter.apply(u, mesh); // warm-up
 
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < repeats; ++i) {
        limiter.apply(u, mesh);
    }
    return msSince(start) / repeats;
}
 
void runLimiterSweep(int repeats) {
    const std::vector<Index> sizes = {128,  256,  512,  1024, 2048,
                                       4096, 8192, 16384, 32768, 65536};
 
    std::printf("=== Part 1b: MinmodLimiter::apply() vs element count ===\n");
    std::printf("%10s %15s\n", "elements", "ms/apply");
    for (Index n : sizes) {
        const double ms = benchmarkLimiter(n, repeats);
        std::printf("%10zu %15.6f\n", static_cast<size_t>(n), ms);
    }
    std::printf("\n");
}
 
// ---------------------------------------------------------------------
// Part 2: full Solver::step() breakdown on a realistic Y-bifurcation.
// ---------------------------------------------------------------------
 
Network makeBifurcationNetwork(Index elementsPerVessel) {
    VesselParameters params;
    params.A0 = 0.126;
    params.beta = 6.06e5;
    params.alpha = 4.0 / 3.0;
    params.frictionKr = 0.0;
    params.nElements = elementsPerVessel;
 
    // Keep element size h ~ 1 regardless of elementsPerVessel (as in
    // makeSingleVesselNetwork above) so the CFL limit on dt doesn't tighten
    // as the benchmark scales up -- otherwise a dt safe for the smallest
    // size becomes unstable for the largest one.
    VesselParameters p1 = params;
    p1.length = static_cast<Real>(elementsPerVessel);
    VesselParameters p2 = params;
    p2.length = static_cast<Real>(elementsPerVessel);
    VesselParameters p3 = params;
    p3.length = static_cast<Real>(elementsPerVessel) * 2.0;
    p3.nElements = elementsPerVessel * 2;
 
    std::vector<Vessel> vessels{Vessel(1, "Omega1", p1), Vessel(2, "Omega2", p2), Vessel(3, "Omega3", p3)};
 
    BoundaryConditionSpec nonReflecting;
    nonReflecting.type = BoundaryConditionType::NonReflecting;
 
    std::vector<Node> nodes{
        Node(1, "inlet", {{1, VesselEnd::Proximal}}, nonReflecting),
        Node(2, "bifurcation",
             {{1, VesselEnd::Distal}, {2, VesselEnd::Proximal}, {3, VesselEnd::Proximal}}),
        Node(3, "outlet2", {{2, VesselEnd::Distal}}, nonReflecting),
        Node(4, "outlet3", {{3, VesselEnd::Distal}}, nonReflecting),
    };
    return Network(FluidProperties{}, std::move(vessels), std::move(nodes));
}
 
struct StepBreakdown {
    double spatialOperatorMs = 0.0; // both RK sub-stage evaluations
    double limiterMs = 0.0;         // both limiter applications
    double boundaryResolveMs = 0.0;
    double combineMs = 0.0; // the two RK stage-combine loops
    double totalStepMs = 0.0;
    Index elementCount = 0;
    Index dofCount = 0;
};
 
StepBreakdown benchmarkStep(Index elementsPerVessel, int repeats) {
    const Network network = makeBifurcationNetwork(elementsPerVessel);
    const Mesh mesh(network);
 
    // A genuine rest state (Q=0, A=A0 everywhere): with all-non-reflecting
    // boundaries and zero friction this is an exact steady fixed point (see
    // the SpatialOperator "uniform state gives zero residual" tests), so it
    // stays numerically stable indefinitely -- unlike a nonzero mean flow,
    // which drifts with nothing anchoring the mean pressure (as found
    // tuning the heartbeat example) and can blow up over many repeated
    // timing iterations. Only wall-clock cost is measured here, not
    // physical behavior, so a stable fixed point is exactly what's wanted.
    State u(mesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = 0.126;
        u.Q[i] = 0.0;
    }
 
    LinearElasticTubeLaw law;
    HllFlux flux;
    MinmodLimiter limiter;
    BoundarySolver resolver(network, mesh, law, network.fluid());
    resolver.solve(u, 0.0, 1e-5);
    Solver solver(mesh, network.fluid(), law, flux, resolver, &limiter);
 
    const Real dt = 1e-5;
 
    // Warm-up.
    for (int i = 0; i < 3; ++i) {
        resolver.solve(u, 0.0, dt);
        solver.step(u, 0.0, dt);
    }
 
    StepBreakdown result;
    result.elementCount = mesh.elementCount();
    result.dofCount = mesh.totalDofs();
 
    // Full, realistic step (resolve + step), as actually used by Simulation::step.
    {
        const auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < repeats; ++i) {
            resolver.solve(u, 0.0, dt);
            solver.step(u, 0.0, dt);
        }
        result.totalStepMs = msSince(start) / repeats;
    }
 
    // Component breakdown, replicating what Solver::step() does internally.
    {
        State k1(mesh.totalDofs()), k2(mesh.totalDofs()), stage(mesh.totalDofs());
        SpatialOperator op(mesh, network.fluid(), law, flux, resolver);
 
        double spatialMs = 0.0, limiterMsAcc = 0.0, combineMsAcc = 0.0, resolveMs = 0.0;
        for (int i = 0; i < repeats; ++i) {
            auto t0 = std::chrono::steady_clock::now();
            resolver.solve(u, 0.0, dt);
            resolveMs += msSince(t0);
 
            t0 = std::chrono::steady_clock::now();
            op.evaluate(u, 0.0, k1);
            spatialMs += msSince(t0);
 
            t0 = std::chrono::steady_clock::now();
            for (Index j = 0; j < u.size(); ++j) {
                stage.A[j] = u.A[j] + dt * k1.A[j];
                stage.Q[j] = u.Q[j] + dt * k1.Q[j];
            }
            combineMsAcc += msSince(t0);
 
            t0 = std::chrono::steady_clock::now();
            limiter.apply(stage, mesh);
            limiterMsAcc += msSince(t0);
 
            t0 = std::chrono::steady_clock::now();
            op.evaluate(stage, dt, k2);
            spatialMs += msSince(t0);
 
            t0 = std::chrono::steady_clock::now();
            for (Index j = 0; j < u.size(); ++j) {
                u.A[j] += 0.5 * dt * (k1.A[j] + k2.A[j]);
                u.Q[j] += 0.5 * dt * (k1.Q[j] + k2.Q[j]);
            }
            combineMsAcc += msSince(t0);
 
            t0 = std::chrono::steady_clock::now();
            limiter.apply(u, mesh);
            limiterMsAcc += msSince(t0);
        }
        result.boundaryResolveMs = resolveMs / repeats;
        result.spatialOperatorMs = spatialMs / repeats;
        result.combineMs = combineMsAcc / repeats;
        result.limiterMs = limiterMsAcc / repeats;
    }
 
    return result;
}
 
void runStepBreakdown(int repeats) {
    std::printf("=== Part 2: full Solver::step() breakdown (Y-bifurcation network) ===\n");
    std::printf("%10s %8s %10s %10s %10s %10s %10s %10s\n", "elements", "dofs", "resolve", "spatialOp",
                "combine", "limiter", "sum", "actualStep");
    for (Index n : {8, 16, 32, 64, 128, 256, 512, 1024, 2048}) {
        const StepBreakdown b = benchmarkStep(n, repeats);
        const double sum = b.boundaryResolveMs + b.spatialOperatorMs + b.combineMs + b.limiterMs;
        std::printf("%10zu %8zu %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
                    static_cast<size_t>(b.elementCount), static_cast<size_t>(b.dofCount),
                    b.boundaryResolveMs, b.spatialOperatorMs, b.combineMs, b.limiterMs, sum,
                    b.totalStepMs);
    }
    std::printf("(all times in ms; 'elements' is total across all 3 vessels)\n");
}
 
} // namespace
 
int main(int argc, char** argv) {
    const int repeats = argc > 1 ? std::atoi(argv[1]) : 200;
    runEvaluateSweep(repeats);
    runLimiterSweep(repeats);
    runStepBreakdown(repeats);
    return 0;
}
 