#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/solver.hpp"
 
using namespace hemo1d;
using namespace hemo1d::dg;
using namespace hemo1d::physics;
using Catch::Approx;
 
namespace {
 
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
 
class ConstantBoundaryStateProvider : public BoundaryStateProvider {
public:
    ConstantBoundaryStateProvider(Real A, Real Q) : A_(A), Q_(Q) {}
    std::pair<Real, Real> exteriorState(Id, VesselEnd, Real) const override { return {A_, Q_}; }
 
private:
    Real A_;
    Real Q_;
};
 
Network makeSingleVesselNetwork() {
    VesselParameters params;
    params.length = 2.0;
    params.A0 = kA0;
    params.beta = kBeta;
    params.frictionKr = 0.0;
    params.nElements = 5;
    params.polynomialOrder = 2;
 
    std::vector<Vessel> vessels{Vessel(1, "v", params)};
    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, bc),
        Node(2, "out", {{1, VesselEnd::Distal}}, bc),
    };
    return Network(FluidProperties{}, std::move(vessels), std::move(nodes));
}
 
} // namespace
 
TEST_CASE("Solver::step leaves an exact steady uniform state unchanged", "[solver]") {
    const Network network = makeSingleVesselNetwork();
    const Mesh mesh(network);
 
    const Real Q0 = 0.05;
    State u(mesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = kA0;
        u.Q[i] = Q0;
    }
    const State initial = u;
 
    LinearElasticTubeLaw law;
    LaxFriedrichsFlux flux;
    ConstantBoundaryStateProvider provider(kA0, Q0);
    MinmodLimiter limiter;
    Solver solver(mesh, FluidProperties{}, law, flux, provider, &limiter);
 
    Real time = 0.0;
    const Real dt = 1e-4;
    for (int step = 0; step < 20; ++step) {
        solver.step(u, time, dt);
        time += dt;
    }
 
    for (Index i = 0; i < u.size(); ++i) {
        CHECK(u.A[i] == Approx(initial.A[i]).margin(1e-8));
        CHECK(u.Q[i] == Approx(initial.Q[i]).margin(1e-8));
    }
}
 
TEST_CASE("Solver::step works without a limiter", "[solver]") {
    const Network network = makeSingleVesselNetwork();
    const Mesh mesh(network);
 
    const Real Q0 = 0.02;
    State u(mesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = kA0;
        u.Q[i] = Q0;
    }
 
    LinearElasticTubeLaw law;
    HllFlux flux;
    ConstantBoundaryStateProvider provider(kA0, Q0);
    Solver solver(mesh, FluidProperties{}, law, flux, provider); // no limiter
 
    solver.step(u, 0.0, 1e-4);
 
    for (Index i = 0; i < u.size(); ++i) {
        CHECK(u.A[i] == Approx(kA0).margin(1e-8));
        CHECK(u.Q[i] == Approx(Q0).margin(1e-8));
    }
}
 