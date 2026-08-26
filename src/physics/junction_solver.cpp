#include "hemo1d/physics/junction_solver.hpp"
#include "hemo1d/core/types.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "hemo1d/core/linear_algebra.hpp"
#include "hemo1d/physics/compatibility.hpp"

namespace hemo1d::physics {

namespace{

// Total pressure formula
Real totalPressure(Real A, Real Q, Real A0, Real beta, Real rho, const TubeLaw& tubeLaw){
    const Real u = Q / A;
    return tubeLaw.pressure(A, A0, beta) + 0.5 * rho * u * u;
}

// Total pressure derivative respect to A
Real totalPressureDA(Real A, Real Q, Real A0, Real beta, Real rho, const TubeLaw& tubeLaw){
    return tubeLaw.pressureDerivative(A, A0, beta) - rho * Q * Q / (A * A * A);
}

// Total pressure derivative respect to Q
Real totalPressureDQ(Real A, Real Q, Real rho) {
    return rho * Q / (A * A);
}

void MassConservationUpdate(const std::vector<Real>& x, const std::vector<Real>& sigma, 
            std::vector<Real>& F, DenseMatrix& J, Index N){
    for (std::size_t i = 0; i < N; ++i) {
        F[0] += sigma[i] * x[2 * i + 1];
        J(0, 2 * i + 1) = sigma[i];
    }
}

void TotalPressureUpdate(const std::vector<Real>& x, const std::vector<Real>& sigma, Real rho, const TubeLaw& tubeLaw,
            std::vector<Real>& F, DenseMatrix& J, const std::vector<JunctionBranch>& branches, Index N){
    
    const Real A0v = x[0];
    const Real Q0v = x[1];
    const Real Ptot0 = totalPressure(A0v, Q0v, branches[0].A0, branches[0].beta, rho, tubeLaw);
    const Real dPtot0dA = totalPressureDA(A0v, Q0v, branches[0].A0, branches[0].beta, rho, tubeLaw);
    const Real dPtot0dQ = totalPressureDQ(A0v, Q0v, rho);
    
    for (std::size_t k=0; k+1<N; ++k){
        const std::size_t i=k+1;
        const Real Ai = x[2*i];
        const Real Qi = x[2*i+1];
        J(i, 0) = dPtot0dA;
        J(i, 1) = dPtot0dQ;

        F[i] = Ptot0 - totalPressure(Ai, Qi, branches[i].A0, branches[i].beta, rho, tubeLaw);
        J(i, 2*i) = -totalPressureDA(Ai, Qi, branches[i].A0, branches[i].beta, rho, tubeLaw);
        J(i, 2*i+1) = -totalPressureDQ(Ai, Qi, rho);
    }
}

void compatibilityUpdate(const std::vector<Real>& x, std::vector<Real>& F, const std::vector<Real>& compatRhs,
        DenseMatrix& J, const std::vector<std::pair<Real, Real>>& lOut, Index N){
        for (std::size_t i=0; i<N; ++i){
            const std::size_t row = N+i;
            F[row] = lOut[i].first * x[2*i] + lOut[i].second * x[2*i+1] - compatRhs[i];
            J(row, 2*i) = lOut[i].first;
            J(row, 2*i+1) = lOut[i].second;
        }
    }

} // namespace



JunctionSolution solveJunction(const std::vector<JunctionBranch>& branches, Real rho,
                    const TubeLaw& tubeLaw, Real dt, Real residualTol, Real incrementTol, int maxIterations){
    
    const std::size_t N = branches.size();
    if (N<2){
        throw std::invalid_argument("solveJunction: at least 2 branches are required");
    }

    std::vector<Real> sigma(N);
    std::vector<std::pair<Real, Real>> lOut(N);
    std::vector<Real> compatRhs(N);

    for (std::size_t i=0; i<N; ++i){
        const JunctionBranch& b = branches[i];
        sigma[i] = (b.end==VesselEnd::Distal) ? 1.0 : -1.0;


        const LeftEigenvectors eig = computeLeftEigenvectors(b.A, b.Q, b.A0, 
            b.beta, b.alpha, rho, tubeLaw);
        
        lOut[i] = (b.end == VesselEnd::Distal) ? eig.lPlus : eig.lMinus;

        const auto [compatA, compatQ] = compatibilityPrediction(b.A, b.Q, b.dAdz, b.dQdz, b.A0,
            b.beta, b.alpha, rho, b.frictionKr, tubeLaw, dt);
        
        compatRhs[i] = lOut[i].first * compatA + lOut[i].second * compatQ;

    }

    const Index n= 2*N;
    std::vector<Real> x(n);
    for (std::size_t i=0; i<N; ++i){
        x[2*i] = branches[i].A;
        x[2*i+1] = branches[i].Q;
    }

    // Newton raphson method
    int iter = 0;
    bool converged = false, stalled = false;
    Real residualNorm = 0.0;
    for (; iter<maxIterations; ++iter){
        std::vector<Real> F(n, 0.0);
        DenseMatrix J(n, n, 0.0);
        MassConservationUpdate(x, sigma, F, J, N);
        TotalPressureUpdate(x, sigma, rho, tubeLaw, F, J, branches, N);
        compatibilityUpdate(x, F, compatRhs, J, lOut, N);

        // residual norm calc:
        Real fNormSq = 0.0;
        for (Real f : F) fNormSq += f*f;
        residualNorm = std::sqrt(fNormSq);
        if (residualNorm < residualTol){
            converged = true;
            break;
        }
        // -1 * F
        std::vector<Real> negF(n);
        for (Index k = 0; k < n; ++k) negF[k] = -F[k];
        const std::vector<Real> delta = solveLinearSystem(J, negF);

        // Condition break
        Real maxRelChange = 0.0;
        for (Index k=0; k<n; ++k){
            x[k] += delta[k];

            maxRelChange = std::max(maxRelChange, std::abs(delta[k]) / 
                                                    (1.0 + std::abs(x[k])));
        }

        if (maxRelChange < incrementTol){
            stalled = true;
            ++iter;
            break;
        }

    }

    JunctionSolution solution;
    solution.A.resize(N);
    solution.Q.resize(N);
    for (std::size_t i=0; i<N; ++i){
        solution.A[i] = x[2*i];
        solution.Q[i] = x[2*i+1];
    }
    solution.iterations = iter;
    solution.converged = converged;
    solution.stalled = stalled;
    solution.residualNorm = residualNorm;
    return solution;
}


} // namespace hemo1d::physics