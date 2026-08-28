#pragma once

#include <vector>

#include "hemo1d/core/vessel.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d::physics{

struct JunctionBranch {
    Real A0 = 0.0;
    Real beta = 0.0;
    Real alpha = 4.0 / 3.0;
    Real frictionKr = 0.0;
    VesselEnd end = VesselEnd::Distal;
    Real A = 0.0;
    Real Q = 0.0;
    Real dAdz = 0.0;
    Real dQdz = 0.0;
};

struct JunctionSolution{
    std::vector<Real> A;
    std::vector<Real> Q;
    int iterations = 0;
    bool converged = false;
    bool stalled = false;
    Real residualNorm = 0.0;
};

JunctionSolution solveJunction(const std::vector<JunctionBranch>& branches, Real rho,
                    const TubeLaw& tubeLaw, Real dt, 
                    Real residualTol = 1e-8, Real incrementTol = 1e-8,int maxIterations = 50);

} // namespace hemo1d::physics