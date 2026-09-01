#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/core/linear_algebra.hpp"
 
using namespace hemo1d;
using Catch::Approx;
 
TEST_CASE("solveLinearSystem solves a small known system", "[linear_algebra]") {
    // [2 1] [x]   [5]      exact solution x=2, y=1
    // [1 3] [y] = [5]
    DenseMatrix A(2, 2);
    A(0, 0) = 2.0;
    A(0, 1) = 1.0;
    A(1, 0) = 1.0;
    A(1, 1) = 3.0;
    const std::vector<Real> x = solveLinearSystem(A, {5.0, 5.0});
 
    REQUIRE(x.size() == 2);
    CHECK(x[0] == Approx(2.0));
    CHECK(x[1] == Approx(1.0));
}
 
TEST_CASE("solveLinearSystem uses partial pivoting when needed", "[linear_algebra]") {
    // Zero on the diagonal forces a row swap.
    DenseMatrix A(2, 2);
    A(0, 0) = 0.0;
    A(0, 1) = 1.0;
    A(1, 0) = 1.0;
    A(1, 1) = 1.0;
    const std::vector<Real> x = solveLinearSystem(A, {3.0, 4.0});
 
    CHECK(x[0] == Approx(1.0));
    CHECK(x[1] == Approx(3.0));
}
 
TEST_CASE("solveLinearSystem rejects a singular matrix", "[linear_algebra]") {
    DenseMatrix A = DenseMatrix::Zero(2, 2);
    A(0, 0) = 1.0;
    A(0, 1) = 2.0;
    A(1, 0) = 2.0;
    A(1, 1) = 4.0; // row 2 = 2 * row 1
    CHECK_THROWS(solveLinearSystem(A, {1.0, 2.0}));
}
 
TEST_CASE("invertDense produces the identity when multiplied by A", "[linear_algebra]") {
    DenseMatrix A(3, 3);
    A(0, 0) = 4;
    A(0, 1) = 3;
    A(0, 2) = 2;
    A(1, 0) = 1;
    A(1, 1) = 5;
    A(1, 2) = 1;
    A(2, 0) = 2;
    A(2, 1) = 1;
    A(2, 2) = 6;
 
    const DenseMatrix inv = invertDense(A);
    for (Index i = 0; i < 3; ++i) {
        for (Index j = 0; j < 3; ++j) {
            Real sum = 0.0;
            for (Index k = 0; k < 3; ++k) sum += A(i, k) * inv(k, j);
            CHECK(sum == Approx(i == j ? 1.0 : 0.0).margin(1e-10));
        }
    }
}