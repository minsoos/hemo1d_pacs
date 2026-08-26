#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/dg/nodal_basis.hpp"
#include "hemo1d/dg/quadrature.hpp"
 
using namespace hemo1d;
using namespace hemo1d::dg;
using Catch::Approx;
 
namespace {
NodalBasis gllBasis(unsigned order) { return NodalBasis(gaussLobattoLegendre(order + 1).points); }
} // namespace
 
TEST_CASE("NodalBasis satisfies the Kronecker delta property at its own nodes", "[nodal_basis]") {
    const NodalBasis basis = gllBasis(4);
    const auto& nodes = basis.nodes();
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const std::vector<Real> l = basis.evaluate(nodes[i]);
        for (std::size_t j = 0; j < nodes.size(); ++j) {
            CHECK(l[j] == Approx(i == j ? 1.0 : 0.0).margin(1e-12));
        }
    }
}
 
TEST_CASE("NodalBasis is a partition of unity everywhere", "[nodal_basis]") {
    const NodalBasis basis = gllBasis(5);
    for (Real x : {-1.0, -0.7, -0.1, 0.0, 0.33, 0.9, 1.0}) {
        const std::vector<Real> l = basis.evaluate(x);
        Real sum = 0.0;
        for (Real li : l) sum += li;
        CHECK(sum == Approx(1.0).margin(1e-12));
    }
}
 
TEST_CASE("NodalBasis differentiation matrix reproduces exact polynomial derivatives",
          "[nodal_basis]") {
    // Order 3 basis: represent u(x) = x^3 - 2x exactly at the nodes and
    // check D*u matches u'(x) = 3x^2 - 2 at every node.
    const NodalBasis basis = gllBasis(3);
    const auto& nodes = basis.nodes();
    const DenseMatrix& D = basis.differentiationMatrix();
 
    std::vector<Real> u(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        u[i] = nodes[i] * nodes[i] * nodes[i] - 2.0 * nodes[i];
    }
 
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        Real du = 0.0;
        for (std::size_t j = 0; j < nodes.size(); ++j) {
            du += D(i, j) * u[j];
        }
        const Real exact = 3.0 * nodes[i] * nodes[i] - 2.0;
        CHECK(du == Approx(exact).margin(1e-10));
    }
}
 
TEST_CASE("NodalBasis evaluateDerivative matches the differentiation matrix at nodes",
          "[nodal_basis]") {
    const NodalBasis basis = gllBasis(4);
    const auto& nodes = basis.nodes();
    const DenseMatrix& D = basis.differentiationMatrix();
 
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const std::vector<Real> d = basis.evaluateDerivative(nodes[i]);
        for (std::size_t j = 0; j < nodes.size(); ++j) {
            CHECK(d[j] == Approx(D(i, j)).margin(1e-10));
        }
    }
}
 
TEST_CASE("NodalBasis rejects duplicate nodes", "[nodal_basis]") {
    CHECK_THROWS(NodalBasis(std::vector<Real>{-1.0, 0.0, 0.0, 1.0}));
}