#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/dg/reference_element.hpp"
 
using namespace hemo1d;
using namespace hemo1d::dg;
using Catch::Approx;
 
TEST_CASE("ReferenceElement has GLL nodes including both endpoints", "[reference_element]") {
    const ReferenceElement ref(3);
    REQUIRE(ref.numNodes() == 4);
    CHECK(ref.nodes().front() == Approx(-1.0));
    CHECK(ref.nodes().back() == Approx(1.0));
}
 
TEST_CASE("ReferenceElement mass matrix is symmetric and its total sums to the interval length",
          "[reference_element]") {
    for (unsigned order = 1; order <= 5; ++order) {
        const ReferenceElement ref(order);
        const DenseMatrix& M = ref.massMatrix();
        const Index n = ref.numNodes();
 
        Real total = 0.0;
        for (Index i = 0; i < n; ++i) {
            for (Index j = 0; j < n; ++j) {
                CHECK(M(i, j) == Approx(M(j, i)).margin(1e-12));
                total += M(i, j);
            }
        }
        // sum_ij M_ij = integral (sum_i l_i)(sum_j l_j) dx = integral 1 dx = 2
        CHECK(total == Approx(2.0).margin(1e-10));
    }
}
 
TEST_CASE("ReferenceElement mass matrix diagonal entries are positive", "[reference_element]") {
    const ReferenceElement ref(4);
    for (Index i = 0; i < ref.numNodes(); ++i) {
        CHECK(ref.massMatrix()(i, i) > 0.0);
    }
}
 
TEST_CASE("ReferenceElement basisAtQuadrature reproduces the partition of unity",
          "[reference_element]") {
    const ReferenceElement ref(3);
    const DenseMatrix& L = ref.basisAtQuadrature();
    for (Index q = 0; q < static_cast<Index>(L.rows()); ++q) {
        Real sum = 0.0;
        for (Index i = 0; i < static_cast<Index>(L.cols()); ++i) sum += L(q, i);
        CHECK(sum == Approx(1.0).margin(1e-12));
    }
}