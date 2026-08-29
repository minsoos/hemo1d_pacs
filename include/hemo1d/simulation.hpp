#pragma once

#include <memory>
#include <string>
#include <vector>

#include "hemo1d/core/network.hpp"
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/output/probe.hpp"
#include "hemo1d/physics/boundary_solver.hpp"
#include "hemo1d/physics/flux.hpp"
#include "hemo1d/physics/slope_limiter.hpp"
#include "hemo1d/physics/solver.hpp"
#include "hemo1d/physics/state.hpp"
#include "hemo1d/physics/tube_law.hpp"

namespace hemo1d {

enum class FluxKind { LaxFriedrichs, Hll };

struct SimulationSettings {
    unsigned defaultPolynomialOrder = 1;
    FluxKind flux = FluxKind::Hll;
    bool useSlopeLimiter = true;
    Real slopeLimiterTolerance = 1e-8;
};

// High level orchestration facade tying together the complete 1D hemodynamics solver pipeline.
class Simulation {
public:
    explicit Simulation(Network network, SimulationSettings settings = {});

    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) = delete;
    Simulation& operator=(Simulation&&) = delete;

    // Set every DOF to the given flow rate and area.
    // If area <= 0.0, it is set to the initial area of the owning vessel.
    void setUniformInitialCondition(Real flowRate = 0.0, Real area = -1.0);

    // Add a probe to a specific vessel of the network, to retrieve the 
    // values of its state at a specific physical point z.
    output::Probe& addProbe(const std::string& name, Id vesselId, Real z);
    const std::vector<output::ProbeSample>& probeSamples(const std::string& name) const;
    void writeProbesCsv(const std::string& directory) const;

    // A real-only reflection of every DOF in the mesh at the current state/time. 
    // Can be used to reconstruct the true discrete solution. 
    std::vector<output::FieldSample> fieldSnapshot() const;

    // Advances the simulation by one step of size dt.
    void step(Real dt);

    // Returns the largest stable value of dt for the current state of the simulation.
    Real cflTimeStep(Real cflNumber = 0.9) const;

    // Advances the simulation until targetTime, recording every probe every `recordEvery` steps.
    // If dt <= 0.0, a dynamic CFL-based dt is used for every step.
    void run(
        Real targetTime, Real dt, int recordEvery = 1, const std::string& vtkDirectory = "",
        int vtkEvery = 0, Real cflNumber = 0.9
    );

    Real time() const noexcept { return time_; }
    const physics::State& state() const noexcept { return state_; }
    const Network& network() const noexcept { return network_; }
    const dg::Mesh& mesh() const noexcept { return mesh_; }

private:
    Network network_;
    dg::Mesh mesh_;
    std::unique_ptr<physics::TubeLaw> tubeLaw_;
    std::unique_ptr<physics::NumericalFlux> flux_;
    std::unique_ptr<physics::MinmodLimiter> limiter_;
    std::unique_ptr<physics::BoundarySolver> boundarySolver_;
    std::unique_ptr<physics::Solver> solver_;
    physics::State state_;
    Real time_ = 0.0;
    output::ProbeRecorder probeRecorder_;
};

}
