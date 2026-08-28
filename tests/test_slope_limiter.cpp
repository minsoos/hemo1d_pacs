#include <cmath>
 
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/physics/slope_limiter.hpp"
 
using namespace hemo1d;
using namespace hemo1d::dg;
using namespace hemo1d::physics;
using Catch::Approx;
 
namespace {
 
Network makeSingleVesselNetwork(Index nElements, unsigned order) {
    VesselParameters params;
    params.length = 4.0;
    params.A0 = 0.126;
    params.beta = 6.06e5;
    params.nElements = nElements;
    params.polynomialOrder = order;
 
    std::vector<Vessel> vessels{Vessel(1, "v", params)};
    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, bc),
        Node(2, "out", {{1, VesselEnd::Distal}}, bc),
    };
    return Network(FluidProperties{}, std::move(vessels), std::move(nodes));
}
 
Real cellAverageOf(const Element& el, const std::vector<Real>& component) {
    const ReferenceElement& ref = el.referenceElement();
    const auto& quad = ref.quadrature();
    const DenseMatrix& L = ref.basisAtQuadrature();
    Real sum = 0.0;
    for (Index q = 0; q < quad.points.size(); ++q) {
        Real value = 0.0;
        for (Index i = 0; i < el.numDofs(); ++i) value += L(q, i) * component[el.dofOffset() + i];
        sum += quad.weights[q] * value;
    }
    return 0.5 * sum;
}
 
} // namespace
 
TEST_CASE("MinmodLimiter leaves a globally linear state unchanged", "[slope_limiter]") {
    const Network network = makeSingleVesselNetwork(5, 2);
    const Mesh mesh(network);
 
    State u(mesh.totalDofs());
    const Real slope = 0.01;
    for (const Element& el : mesh.elements()) {
        for (Index i = 0; i < el.numDofs(); ++i) {
            u.A[el.dofOffset() + i] = 0.126 + slope * el.physicalNodes()[i];
            u.Q[el.dofOffset() + i] = 0.05;
        }
    }
    const State original = u;
 
    MinmodLimiter limiter;
    limiter.apply(u, mesh);
 
    for (Index i = 0; i < u.size(); ++i) {
        CHECK(u.A[i] == Approx(original.A[i]).margin(1e-10));
        CHECK(u.Q[i] == Approx(original.Q[i]).margin(1e-10));
    }
}
 
TEST_CASE("MinmodLimiter flattens a spurious spike while preserving the cell average",
          "[slope_limiter]") {
    const Network network = makeSingleVesselNetwork(5, 2);
    const Mesh mesh(network);
 
    State u(mesh.totalDofs());
    for (Index i = 0; i < u.size(); ++i) {
        u.A[i] = 0.126;
        u.Q[i] = 0.05;
    }
 
    // Inject an artificial oscillation into the middle element's interior DOF.
    const Element& spikeElement = mesh.elements()[2];
    u.A[spikeElement.dofOffset() + 1] += 0.05; // interior node of the order-2 element
    const Real averageBeforeLimiting = cellAverageOf(spikeElement, u.A);
 
    MinmodLimiter limiter;
    limiter.apply(u, mesh);
 
    // The spike should have been detected and the element replaced by a
    // linear reconstruction, but the cell average must be preserved exactly.
    const Real averageAfter = cellAverageOf(spikeElement, u.A);
    CHECK(averageAfter == Approx(averageBeforeLimiting).margin(1e-10));
    // And the limiter should indeed have changed something (the spike is real).
    CHECK(std::abs(u.A[spikeElement.dofOffset() + 1] - (0.126 + 0.05)) > 1e-6);
 
    // Neighboring, still-uniform elements should be untouched.
    const Element& farElement = mesh.elements()[0];
    for (Index i = 0; i < farElement.numDofs(); ++i) {
        CHECK(u.A[farElement.dofOffset() + i] == Approx(0.126).margin(1e-10));
    }
}
 