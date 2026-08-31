#pragma once

#include <vector>

#include "hemo1d/core/vessel.hpp"
#include "hemo1d/physics/section_state.hpp"

namespace hemo1d::physics{

class BloodFlowModel;

struct JunctionBranch {
    VesselParameters params{};
    VesselEnd end = VesselEnd::Distal;
    SectionState trace{};
    SectionGradient grad{};
};

struct JunctionSolution{
    std::vector<Real> A;
    std::vector<Real> Q;
    int iterations = 0;
    bool converged = false;
    bool stalled = false;
    Real residualNorm = 0.0;
};

JunctionSolution solveJunction(
    const std::vector<JunctionBranch>& branches, const BloodFlowModel& model,
    Real dt, Real tolerance = 1e-8, int maxIterations = 50
);

} // namespace hemo1d::physics