#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/blood_flow_model.hpp"
#include "hemo1d/physics/cfl.hpp"

using namespace hemo1d;
using namespace hemo1d::dg;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {

constexpr Real kA0 = 0.126;
constexpr Real kBeta = 6.06e5;
constexpr Real kAlpha = 4.0 / 3.0;

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

TEST_CASE("cflTimeStep matches h/((2p+1)*c) for a uniform, at-rest state", "[cfl]") {
    const Index nElements = 6;
    const unsigned order = 2;
    const Network network = makeSingleVesselNetwork(nElements);
    DgSettings settings;
    settings.defaultPolynomialOrder = order;
    const Mesh mesh(network, settings);

    State u(mesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = kA0;
        u.Q[i] = 0.0; // at rest: lambda = +/- c, no alpha*u contribution
    }

    LinearElasticTubeLaw law;
    const FluidProperties fluid{};
    const BloodFlowModel model(law, fluid);
    VesselParameters wall;
    wall.A0 = kA0;
    wall.beta = kBeta;
    const Real c = law.waveSpeed(kA0, wall, fluid.density);
    const Real h = 2.0 / static_cast<Real>(nElements); // vessel length 2.0, uniform mesh
    const Real cflNumber = 0.7;
    const Real expected = cflNumber * h / ((2.0 * order + 1.0) * c);

    CHECK(cflTimeStep(u, mesh, model, cflNumber) == Approx(expected).epsilon(1e-9));
}

TEST_CASE("cflTimeStep scales inversely with cflNumber and is smaller for higher order", "[cfl]") {
    const Network network = makeSingleVesselNetwork(6);
    LinearElasticTubeLaw law;
    const FluidProperties fluid{};
    const BloodFlowModel model(law, fluid);

    DgSettings orderOneSettings;
    orderOneSettings.defaultPolynomialOrder = 1;
    const Mesh meshOrderOne(network, orderOneSettings);
    DgSettings orderThreeSettings;
    orderThreeSettings.defaultPolynomialOrder = 3;
    const Mesh meshOrderThree(network, orderThreeSettings);

    State u1(meshOrderOne.totalDofs());
    for (Index i = 0; i < u1.size(); ++i) {
        u1.A[i] = kA0;
        u1.Q[i] = 0.02;
    }
    State u3(meshOrderThree.totalDofs());
    for (Index i = 0; i < u3.size(); ++i) {
        u3.A[i] = kA0;
        u3.Q[i] = 0.02;
    }

    const Real dtOrderOne = cflTimeStep(u1, meshOrderOne, model, 0.9);
    const Real dtOrderThree = cflTimeStep(u3, meshOrderThree, model, 0.9);
    CHECK(dtOrderThree < dtOrderOne);

    const Real dtHalfCfl = cflTimeStep(u1, meshOrderOne, model, 0.45);
    CHECK(dtHalfCfl == Approx(0.5 * dtOrderOne).epsilon(1e-9));
}
