#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/flux.hpp"
#include "hemo1d/physics/spatial_operator.hpp"
 
using namespace hemo1d;
using namespace hemo1d::dg;
using namespace hemo1d::physics;
using Catch::Approx;
 
namespace {
 
constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kAlpha = 4.0 / 3.0;
 
class ConstantBoundaryStateProvider : public BoundaryStateProvider {
public:
    ConstantBoundaryStateProvider(Real A, Real Q) : A_(A), Q_(Q) {}
    std::pair<Real, Real> exteriorState(Id, VesselEnd, Real) const override { return {A_, Q_}; }
 
private:
    Real A_;
    Real Q_;
};
 
Network makeSingleVesselNetwork(Index nElements) {
    VesselParameters params;
    params.length = 2.0;
    params.A0 = kA0;
    params.beta = kBeta;
    params.alpha = kAlpha;
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
 
} // namespace
 
TEMPLATE_TEST_CASE("SpatialOperator gives (near) zero residual for a uniform state", "[spatial_operator]",
                    LaxFriedrichsFlux, HllFlux) {
    const Network network = makeSingleVesselNetwork(6);
    DgSettings settings;
    settings.defaultPolynomialOrder = 3;
    const Mesh mesh(network, settings);
 
    const Real Q0 = 0.05;
    State u(mesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = kA0;
        u.Q[i] = Q0;
    }
 
    LinearElasticTubeLaw law;
    TestType flux;
    ConstantBoundaryStateProvider provider(kA0, Q0);
    SpatialOperator op(mesh, FluidProperties{}, law, flux, provider);
 
    State dudt(mesh.totalDofs());
    op.evaluate(u, 0.0, dudt);
 
    for (Index i = 0; i < dudt.size(); ++i) {
        CHECK(dudt.A[i] == Approx(0.0).margin(1e-9));
        CHECK(dudt.Q[i] == Approx(0.0).margin(1e-9));
    }
}
 
TEST_CASE("SpatialOperator friction source is uniform for a uniform state", "[spatial_operator]") {
    const Real Q0 = 0.05;
    const Real Kr = 2.0;
    VesselParameters params;
    params.length = 2.0;
    params.A0 = kA0;
    params.beta = kBeta;
    params.alpha = kAlpha;
    params.frictionKr = Kr;
    params.nElements = 4;
    std::vector<Vessel> vessels{Vessel(1, "v", params)};
    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, bc),
        Node(2, "out", {{1, VesselEnd::Distal}}, bc),
    };
    const Network frictionNetwork(FluidProperties{}, std::move(vessels), std::move(nodes));
    const Mesh frictionMesh(frictionNetwork);
 
    State u(frictionMesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = kA0;
        u.Q[i] = Q0;
    }
 
    LinearElasticTubeLaw law;
    LaxFriedrichsFlux flux;
    ConstantBoundaryStateProvider provider(kA0, Q0);
    SpatialOperator op(frictionMesh, FluidProperties{}, law, flux, provider);
 
    State dudt(frictionMesh.totalDofs());
    op.evaluate(u, 0.0, dudt);
 
    const Real expected = -Kr * Q0 / kA0;
    for (Index i = 0; i < dudt.size(); ++i) {
        CHECK(dudt.A[i] == Approx(0.0).margin(1e-9));
        CHECK(dudt.Q[i] == Approx(expected).margin(1e-9));
    }
}