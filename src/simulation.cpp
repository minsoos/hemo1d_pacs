#include "hemo1d/simulation.hpp"

#include <algorithm>

#include "hemo1d/output/vtk_writer.hpp"
#include "hemo1d/physics/cfl.hpp"

namespace hemo1d {

namespace {

std::unique_ptr<physics::NumericalFlux> makeFlux(FluxKind kind) {
    if (kind == FluxKind::Hll) return std::make_unique<physics::HllFlux>();
    return std::make_unique<physics::LaxFriedrichsFlux>();
}

} // namespace

Simulation::Simulation(
    Network network, SimulationSettings settings
)
    : network_(std::move(network)),
    mesh_(network_, dg::DgSettings{settings.defaultPolynomialOrder}),
    tubeLaw_(std::make_unique<physics::LinearElasticTubeLaw>()),
    model_(*tubeLaw_, network_.fluid()),
    flux_(makeFlux(settings.flux)),
    limiter_(
        settings.useSlopeLimiter
            ? std::make_unique<physics::MinmodLimiter>(settings.slopeLimiterTolerance)
            : nullptr
    ),
    boundarySolver_(
        std::make_unique<physics::BoundarySolver>(
            network_, mesh_, *tubeLaw_, network_.fluid()
        )
    ),
    solver_(
        std::make_unique<physics::Solver>(
            mesh_, model_, *flux_,
            *boundarySolver_, limiter_.get()
        )
    ),
    state_(mesh_.totalDofs())
{
    setUniformInitialCondition();
}

void Simulation::setUniformInitialCondition(
    Real flowRate, Real area
) {
    for (const dg::Element& el : mesh_.elements()) {
        const Real A = area > 0.0 ? area : el.vesselParameters().A0;
        for (Index i = 0; i < el.numDofs(); ++i) {
            state_.A[el.dofOffset() + i] = A;
            state_.Q[el.dofOffset() + i] = flowRate;
        }
    }
}

output::Probe& Simulation::addProbe(
    const std::string& name, Id vesselId, Real z
) {
    return probeRecorder_.addProbe(name, vesselId, z, mesh_);
}

const std::vector<output::ProbeSample>& Simulation::probeSamples(
    const std::string& name
) const {
    return probeRecorder_.samples(name);
}

void Simulation::writeProbesCsv(
    const std::string& directory
) const {
    probeRecorder_.writeCsv(directory);
}

std::vector<output::FieldSample> Simulation::fieldSnapshot() const {
    std::vector<output::FieldSample> samples;
    samples.reserve(state_.size());
    for (const dg::Element& el : mesh_.elements()) {
        const std::vector<Real>& nodes = el.physicalNodes();
        const Index offset = el.dofOffset();
        for (Index i = 0; i < el.numDofs(); ++i) {
            samples.push_back({
                el.flatIndex(), el.vesselId(), nodes[i], 
                state_.A[offset + i], state_.Q[offset + i]
            });
        }
    }
    return samples;
}

void Simulation::step(Real dt) {
    boundarySolver_->solve(state_, time_, dt);
    solver_->step(state_, time_, dt);
    time_ += dt;
}

Real Simulation::cflTimeStep(Real cflNumber) const {
    return physics::cflTimeStep(state_, mesh_, model_, cflNumber);
}

void Simulation::run(
    Real targetTime, Real dt, int recordEvery, const std::string& vtkDirectory,
    int vtkEvery, Real cflNumber
) {
    std::unique_ptr<output::VtkSeriesWriter> vtkWriter;
    if (!vtkDirectory.empty()) {
        vtkWriter = std::make_unique<output::VtkSeriesWriter>(vtkDirectory, "hemo1d");
    }

    const bool dynamicDt = dt <= 0.0;
    constexpr Real kEndTolerance = 1e-12;
    int stepIndex = 0;

    while (time_ < targetTime - kEndTolerance) {
        const Real requestedDt = dynamicDt ? cflTimeStep(cflNumber) : dt;
        const Real thisDt = std::min(requestedDt, targetTime - time_);
        step(thisDt);
        ++stepIndex;

        if (recordEvery > 0 && stepIndex % recordEvery == 0) {
            probeRecorder_.record(state_, time_, *tubeLaw_);
        }

        if (vtkWriter && vtkEvery > 0 && stepIndex % vtkEvery == 0) {
            vtkWriter->write(mesh_, state_, *tubeLaw_, network_.fluid().density, time_);
        }
    }
}

} // namespace hemo1d
