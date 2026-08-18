#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/core/network.hpp"
 
using namespace hemo1d;
 
namespace {
 
Vessel makeVessel(Id id) {
    VesselParameters p;
    p.length = 1.0;
    p.A0 = 0.1;
    p.beta = 1.0;
    return Vessel(id, "v" + std::to_string(id), p);
}
 
BoundaryConditionSpec nonReflecting() {
    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;
    return bc;
}

} // namespace

TEST_CASE("Network: Simplest network") {
    std::vector<Vessel> vessels{makeVessel(1)};
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, nonReflecting()),
        Node(2, "out", {{1, VesselEnd::Distal}}, nonReflecting()),
    };
 
    const Network net(FluidProperties{}, std::move(vessels), std::move(nodes));
    CHECK(net.vesselCount() == 1);
    CHECK(net.nodeCount() == 2);
    CHECK(net.nodeIdAt(1, VesselEnd::Proximal) == 1);
    CHECK(net.nodeIdAt(1, VesselEnd::Distal) == 2);
}

TEST_CASE("Network: Duplicated vessels") {
    std::vector<Vessel> vessels{makeVessel(1), makeVessel(1)};
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, nonReflecting()),
        Node(2, "out", {{1, VesselEnd::Distal}}, nonReflecting()),
    };
 
    CHECK_THROWS(Network(FluidProperties{}, std::move(vessels), std::move(nodes)));
}

TEST_CASE("Network: Duplicated nodes") {
    std::vector<Vessel> vessels{makeVessel(1)};
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, nonReflecting()),
        Node(1, "out", {{1, VesselEnd::Distal}}, nonReflecting()),
    };
 
    CHECK_THROWS(Network(FluidProperties{}, std::move(vessels), std::move(nodes)));
}

TEST_CASE("Network: Duplicated node type") {
    std::vector<Vessel> vessels{makeVessel(1)};
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, nonReflecting()),
        Node(2, "out", {{1, VesselEnd::Proximal}}, nonReflecting()),
    };
 
    CHECK_THROWS(Network(FluidProperties{}, std::move(vessels), std::move(nodes)));
}

TEST_CASE("Network: 2 nodes assigning to 1 space") {
    std::vector<Vessel> vessels{makeVessel(1)};
    std::vector<Node> nodes{
        Node(1, "in1", {{1, VesselEnd::Proximal}}, nonReflecting()),
        Node(2, "in2", {{1, VesselEnd::Proximal}}, nonReflecting()),
        Node(3, "out", {{1, VesselEnd::Distal}}, nonReflecting()),
    };
 
    CHECK_THROWS(Network(FluidProperties{}, std::move(vessels), std::move(nodes)));
}

TEST_CASE("Network: Juntion creation") {
    std::vector<Vessel> vessels{makeVessel(1), makeVessel(2), makeVessel(3)};
    std::vector<Node> nodes{
        Node(1, "in1", {{1, VesselEnd::Proximal}}, nonReflecting()),
        Node(2, "in2", {{2, VesselEnd::Proximal}}, nonReflecting()),
        Node(3, "junction",
             {{1, VesselEnd::Distal}, {2, VesselEnd::Distal}, {3, VesselEnd::Proximal}},
             std::nullopt, {0.5, 1.0}),
        Node(4, "out", {{3, VesselEnd::Distal}}, nonReflecting()),
    };
    const Network net(FluidProperties{}, std::move(vessels), std::move(nodes));
    CHECK(net.node(3).kind() == NodeKind::Junction);
}