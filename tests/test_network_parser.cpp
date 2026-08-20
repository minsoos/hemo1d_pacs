#include <filesystem>
 
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "hemo1d/io/network_parser.hpp"

using namespace hemo1d;
using Catch::Approx;

namespace{
std::filesystem::path examplePath(const std::string& name){
    return std::filesystem::path(HEMO1D_EXAMPLES_DIR) / name;
}
} // namespace

TEST_CASE("Network Parser: loadNetwork: Simple Bifurcation example", "[network_parser]"){
    const Network net = io::loadNetwork(examplePath("simple_bifurcation.json"));

CHECK(net.vesselCount() == 3);
    CHECK(net.nodeCount() == 4);
 
    CHECK(net.fluid().density == Approx(1.055));
    CHECK(net.fluid().viscosity == Approx(0.045));
 
    const Vessel& omega1 = net.vessel(1);
    CHECK(omega1.name() == "Omega1");
    CHECK(omega1.parameters().length == Approx(1.0));
    CHECK(omega1.parameters().A0 == Approx(0.126));
    CHECK(omega1.parameters().beta == Approx(606060.0));
    CHECK(omega1.parameters().nElements == 32);
 
    const Vessel& omega3 = net.vessel(3);
    CHECK(omega3.parameters().length == Approx(2.0));
    CHECK(omega3.parameters().nElements == 64);
 
    // Topology: vessel 1 distal end and vessel 2/3 proximal ends meet at the
    // bifurcation node (id 2).
    CHECK(net.nodeIdAt(1, VesselEnd::Distal) == 2);
    CHECK(net.nodeIdAt(2, VesselEnd::Proximal) == 2);
    CHECK(net.nodeIdAt(3, VesselEnd::Proximal) == 2);
 
    const Node& bifurcation = net.node(2);
    CHECK(bifurcation.kind() == NodeKind::Junction);
    CHECK(bifurcation.connections().size() == 3);
    REQUIRE(bifurcation.bifurcationAnglesRad().size() == 2);
    CHECK(bifurcation.bifurcationAnglesRad()[0] == Approx(0.7853981634));
 
    const Node& inlet = net.node(1);
    CHECK(inlet.kind() == NodeKind::Terminal);
    REQUIRE(inlet.boundaryCondition().has_value());
    CHECK(inlet.boundaryCondition()->type == BoundaryConditionType::Prescribed);
    CHECK(inlet.boundaryCondition()->quantity == PrescribedQuantity::Area);
    CHECK(std::filesystem::exists(inlet.boundaryCondition()->csvFile));
 
    const Node& outlet3 = net.node(4);
    REQUIRE(outlet3.boundaryCondition().has_value());
    CHECK(outlet3.boundaryCondition()->type == BoundaryConditionType::NonReflecting);

}

 
TEST_CASE("Network Parser: loadNetwork: Throws an error on a missing file", "[network_parser]") {
    CHECK_THROWS(io::loadNetwork("does/not/exist.json"));
}