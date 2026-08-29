#include <filesystem>
#include <fstream>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/output/probe.hpp"

using namespace hemo1d;
using namespace hemo1d::dg;
using namespace hemo1d::output;
using namespace hemo1d::physics;
using Catch::Approx;

namespace {

Network makeSingleVesselNetwork() {
VesselParameters params;
params.length = 4.0;
params.A0 = 0.126;
params.beta = 6.06e5;
params.nElements = 4;
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

TEST_CASE("Probe: Single vessel with a prescribed inlet and stays stable",
          "[probe]") {
    const Network network = makeSingleVesselNetwork();
    const Mesh mesh(network);
 
    State state(mesh.totalDofs());
    for (Index i = 0; i < state.size(); ++i) {
        state.A[i] = 0.126;
        state.Q[i] = 0.05;
    }
 
    LinearElasticTubeLaw law;
    for (Real z : {0.0, 0.37, 1.0, 2.5, 3.9999, 4.0}) {
        Probe probe("p", 1, z, mesh);
        const ProbeSample s = probe.sample(state, 1.23, law);
        CHECK(s.time == Approx(1.23));
        CHECK(s.A == Approx(0.126).margin(1e-10));
        CHECK(s.Q == Approx(0.05).margin(1e-10));
        CHECK(s.velocity == Approx(0.05 / 0.126).margin(1e-8));
        CHECK(s.pressure == Approx(0.0).margin(1e-8));
    }

}

} // namespace