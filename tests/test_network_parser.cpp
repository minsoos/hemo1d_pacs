#include <filesystem>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "hemo1d/io/network_parser.hpp"
#include "hemo1d/couplings/windkessel_coupling.hpp"

using namespace hemo1d;
using Catch::Approx;

namespace {
std::filesystem::path examplePath(const std::string& name) {
    return std::filesystem::path(HEMO1D_EXAMPLES_DIR) / name;
}
} // namespace

TEST_CASE("loadNetwork parses the simple bifurcation example", "[network_parser]") {
    const Network net = io::loadNetwork(examplePath("simple_bifurcation.json"));

    CHECK(net.vesselCount() == 3);
    CHECK(net.nodeCount() == 4);

    CHECK(net.fluid().density == Approx(1.05));
    CHECK(net.fluid().viscosity == Approx(0.035));

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

TEST_CASE("loadNetwork parses an external (windkessel) terminal", "[network_parser]") {
    const Network net = io::loadNetwork(examplePath("single_vessel_windkessel.json"));

    REQUIRE(net.nodeCount() == 2);
    const Node& outlet = net.node(2);
    REQUIRE(outlet.kind() == NodeKind::Terminal);
    REQUIRE(outlet.boundaryCondition().has_value());
    CHECK(outlet.boundaryCondition()->type == BoundaryConditionType::External);
    CHECK(outlet.boundaryCondition()->modelName == "windkessel");

    // Core keeps the params as an opaque JSON string; the model parses it.
    const couplings::WindkesselParameters wk =
        couplings::parseWindkesselParams(outlet.boundaryCondition()->modelParams);
    CHECK(wk.R1 < 0.0); // matched-impedance sentinel
    REQUIRE(wk.compartments.size() == 2);
    CHECK(wk.compartments[0].resistance == Approx(5.0e4));
    CHECK(wk.compartments[0].compliance == Approx(8.0e-7));
    CHECK(wk.compartments[1].resistance == Approx(3.0e4));
    CHECK(wk.pOut == Approx(0.0));
    CHECK(wk.pOutCsv.empty());
    CHECK(wk.subSteps == 2);
}

TEST_CASE("loadNetwork throws a helpful error on a missing file", "[network_parser]") {
    CHECK_THROWS(io::loadNetwork("does/not/exist.json"));
}
