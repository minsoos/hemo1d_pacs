#include "hemo1d/physics/junction_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "hemo1d/core/linear_algebra.hpp"
#include "hemo1d/physics/blood_flow_model.hpp"

namespace hemo1d::physics {

namespace {

// Derivatives of BloodFlowModel::totalPressure w.r.t. the branch state,
// needed for the Newton Jacobian of the pressure-continuity rows.
Real totalPressureDA(SectionState u, const VesselParameters& p, const BloodFlowModel& model) {
    return model.tubeLaw().pressureDerivative(u.A, p) -
           model.density() * u.Q * u.Q / (u.A * u.A * u.A);
}

Real totalPressureDQ(SectionState u, Real rho) { return rho * u.Q / (u.A * u.A); }

} // namespace

JunctionSolution solveJunction(const std::vector<JunctionBranch>& branches,
                                const BloodFlowModel& model, Real dt, Real tolerance,
                                int maxIterations) {
    const std::size_t N = branches.size();
    if (N < 2) {
        throw std::invalid_argument("solveJunction: at least 2 branches are required");
    }
    const Real rho = model.density();

    // Frozen (explicit, evaluated at the given branch data) quantities: the
    // sign convention for mass conservation, the outgoing left eigenvector,
    // and the compatibility condition's right-hand side.
    std::vector<Real> sigma(N);
    std::vector<std::pair<Real, Real>> lOut(N);
    std::vector<Real> compatRhs(N);

    for (std::size_t i = 0; i < N; ++i) {
        const JunctionBranch& b = branches[i];
        sigma[i] = (b.end == VesselEnd::Distal) ? 1.0 : -1.0;

        const LeftEigenvectors eig = model.leftEigenvectors(b.trace, b.params);
        lOut[i] = (b.end == VesselEnd::Proximal) ? eig.lMinus : eig.lPlus;

        const SectionState cc = model.compatibilityPrediction(b.trace, b.grad, b.params, dt);
        compatRhs[i] = lOut[i].first * cc.A + lOut[i].second * cc.Q;
    }

    const Index n = 2 * N;
    std::vector<Real> x(n);
    for (std::size_t i = 0; i < N; ++i) {
        x[2 * i] = branches[i].trace.A;
        x[2 * i + 1] = branches[i].trace.Q;
    }

    int iter = 0;
    Real residualNorm = 0.0;
    for (; iter < maxIterations; ++iter) {
        std::vector<Real> F(n, 0.0);
        DenseMatrix J = DenseMatrix::Zero(n, n);

        // Conservation of mass (row 0), linear.
        for (std::size_t i = 0; i < N; ++i) {
            F[0] += sigma[i] * x[2 * i + 1];
            J(0, 2 * i + 1) = sigma[i];
        }

        // Continuity of total pressure, branch 0 vs every other branch
        // (rows 1..N-1), the only nonlinear equations.
        const SectionState u0{x[0], x[1]};
        const Real Ptot0 = model.totalPressure(u0, branches[0].params);
        const Real dPtot0dA = totalPressureDA(u0, branches[0].params, model);
        const Real dPtot0dQ = totalPressureDQ(u0, rho);
        for (std::size_t k = 0; k + 1 < N; ++k) {
            const std::size_t i = k + 1;
            const SectionState ui{x[2 * i], x[2 * i + 1]};
            const std::size_t row = 1 + k;
            F[row] = Ptot0 - model.totalPressure(ui, branches[i].params);
            J(row, 0) = dPtot0dA;
            J(row, 1) = dPtot0dQ;
            J(row, 2 * i) = -totalPressureDA(ui, branches[i].params, model);
            J(row, 2 * i + 1) = -totalPressureDQ(ui, rho);
        }

        // Compatibility conditions (rows N..2N-1), linear.
        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t row = N + i;
            F[row] = lOut[i].first * x[2 * i] + lOut[i].second * x[2 * i + 1] - compatRhs[i];
            J(row, 2 * i) = lOut[i].first;
            J(row, 2 * i + 1) = lOut[i].second;
        }

        Real fNormSq = 0.0;
        for (Real f : F) fNormSq += f * f;
        residualNorm = std::sqrt(fNormSq);

        std::vector<Real> negF(n);
        for (Index k = 0; k < n; ++k) negF[k] = -F[k];
        const std::vector<Real> delta = solveLinearSystem(J, negF);

        Real maxRelChange = 0.0;
        for (Index k = 0; k < n; ++k) {
            x[k] += delta[k];
            maxRelChange = std::max(maxRelChange, std::abs(delta[k]) / (1.0 + std::abs(x[k])));
        }

        if (maxRelChange < tolerance) {
            ++iter;
            break;
        }
    }

    JunctionSolution solution;
    solution.A.resize(N);
    solution.Q.resize(N);
    for (std::size_t i = 0; i < N; ++i) {
        solution.A[i] = x[2 * i];
        solution.Q[i] = x[2 * i + 1];
    }
    solution.iterations = iter;
    solution.residualNorm = residualNorm;
    return solution;
}

} // namespace hemo1d::physics
