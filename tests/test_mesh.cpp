#include <filesystem>
 
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/io/network_parser.hpp"
 
using namespace hemo1d;
using namespace hemo1d::dg;
using Catch::Approx;
 
namespace {
Network loadSimpleBifurcation() {
    return io::loadNetwork(std::filesystem::path(HEMO1D_EXAMPLES_DIR) / "simple_bifurcation.json");
}
} // namespace
 
TEST_CASE("Mesh flattens all vessels into one contiguous element array", "[mesh]") {
    const Network network = loadSimpleBifurcation();
    DgSettings settings;
    settings.defaultPolynomialOrder = 2;
    const Mesh mesh(network, settings);
 
    // n_elements: Omega1=32, Omega2=32, Omega3=64 (see simple_bifurcation.json)
    CHECK(mesh.elementCount() == 32 + 32 + 64);
    CHECK(mesh.totalDofs() == mesh.elementCount() * 3); // order 2 -> 3 nodes per element
 
    const auto [begin1, end1] = mesh.vesselElementRange(1);
    const auto [begin2, end2] = mesh.vesselElementRange(2);
    const auto [begin3, end3] = mesh.vesselElementRange(3);
    CHECK(begin1 == 0);
    CHECK(end1 == 32);
    CHECK(begin2 == 32);
    CHECK(end2 == 64);
    CHECK(begin3 == 64);
    CHECK(end3 == 128);
}
 
TEST_CASE("Mesh elements carry correct neighbor topology within a vessel", "[mesh]") {
    const Network network = loadSimpleBifurcation();
    const Mesh mesh(network);
    const auto [begin, end] = mesh.vesselElementRange(1);
 
    const Element& first = mesh.elements()[begin];
    CHECK(first.isFirstInVessel());
    CHECK_FALSE(first.isLastInVessel());
    REQUIRE(first.rightNeighbor().has_value());
    CHECK(*first.rightNeighbor() == begin + 1);
 
    const Element& last = mesh.elements()[end - 1];
    CHECK(last.isLastInVessel());
    CHECK_FALSE(last.isFirstInVessel());
    REQUIRE(last.leftNeighbor().has_value());
    CHECK(*last.leftNeighbor() == end - 2);
 
    const Element& mid = mesh.elements()[begin + 5];
    REQUIRE(mid.leftNeighbor().has_value());
    REQUIRE(mid.rightNeighbor().has_value());
    CHECK(*mid.leftNeighbor() == begin + 4);
    CHECK(*mid.rightNeighbor() == begin + 6);
}
 
TEST_CASE("Mesh element geometry matches a uniform partition of the vessel", "[mesh]") {
    const Network network = loadSimpleBifurcation();
    const Mesh mesh(network);
    const auto [begin, end] = mesh.vesselElementRange(3); // Omega3: L=2, n_elements=64
 
    const Real h = 2.0 / 64.0;
    for (Index e = begin; e < end; ++e) {
        const Element& el = mesh.elements()[e];
        const Index local = e - begin;
        CHECK(el.zLeft() == Approx(local * h));
        CHECK(el.zRight() == Approx((local + 1) * h));
        CHECK(el.jacobian() == Approx(h / 2.0));
 
        // GLL nodes include the endpoints, so the physical endpoints are
        // exactly the first/last interpolation nodes.
        CHECK(el.physicalNodes().front() == Approx(el.zLeft()));
        CHECK(el.physicalNodes().back() == Approx(el.zRight()));
    }
}
 
TEST_CASE("Mesh assigns contiguous, non-overlapping DOF offsets", "[mesh]") {
    const Network network = loadSimpleBifurcation();
    const Mesh mesh(network);
 
    Index expectedOffset = 0;
    for (const Element& el : mesh.elements()) {
        CHECK(el.dofOffset() == expectedOffset);
        expectedOffset += el.numDofs();
    }
    CHECK(expectedOffset == mesh.totalDofs());
}
 
TEST_CASE("Mesh caches one ReferenceElement per distinct polynomial order", "[mesh]") {
    const Network network = loadSimpleBifurcation();
    const Mesh mesh(network); // all vessels use the default order (1, no override in the example)
 
    CHECK(mesh.referenceElement(1).numNodes() == 2);
    CHECK_THROWS(mesh.referenceElement(7));
}
 
TEST_CASE("Mesh elements cache the owning vessel's physical parameters", "[mesh]") {
    const Network network = loadSimpleBifurcation();
    const Mesh mesh(network);
    const auto [begin, end] = mesh.vesselElementRange(1);
    for (Index e = begin; e < end; ++e) {
        CHECK(mesh.elements()[e].vesselParameters().A0 == Approx(0.126));
        CHECK(mesh.elements()[e].vesselId() == 1);
    }
}
 