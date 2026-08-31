#include "hemo1d/output/vtk_writer.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "hemo1d/physics/blood_flow_model.hpp"

namespace hemo1d::output {

void writeVtk(
    const std::filesystem::path& path, const dg::Mesh& mesh, 
    const physics::State& state, const physics::BloodFlowModel& model, 
    Real time
) {
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("writeVtk: cannot open file " + path.string());
    }
    file << std::setprecision(10);

    const std::vector<dg::Element>& elements = mesh.elements();
    const Index totalPoints = mesh.totalDofs();

    file << "# vtk DataFile Version 3.0\n";
    file << "Hemo1D field data at t=" << time << "\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    constexpr Real kVesselSpacing = 1.0; // visual vertical separation

    file << "POINTS " << totalPoints << " double\n";
    for (const dg::Element& el : elements) {
        const Real y = static_cast<Real>(el.vesselId()) * kVesselSpacing;
        for (Real z : el.physicalNodes()) {
            file << z << ' ' << y << " 0\n";
        }
    }


    Index lineListSize = 0;
    for (const dg::Element& el : elements) lineListSize += el.numDofs() + 1;
    file << "LINES " << elements.size() << ' ' << lineListSize << "\n";
    for (const dg::Element& el : elements) {
        file << el.numDofs();
        for (Index i = 0; i < el.numDofs(); ++i) file << ' ' << (el.dofOffset() + i);
        file << "\n";
    }

    file << "POINT_DATA " << totalPoints << "\n";

    file << "SCALARS Area double 1\nLOOKUP_TABLE default\n";
    for (Real a : state.A) file << a << "\n";

    file << "SCALARS FlowRate double 1\nLOOKUP_TABLE default\n";
    for (Real q : state.Q) file << q << "\n";

    file << "SCALARS Velocity double 1\nLOOKUP_TABLE default\n";
    for (Index i = 0; i < totalPoints; ++i) file << (state.Q[i] / state.A[i]) << "\n";

    file << "SCALARS Pressure double 1\nLOOKUP_TABLE default\n";
    for (const dg::Element& el : elements) {
        const VesselParameters& p = el.vesselParameters();
        for (Index i = 0; i < el.numDofs(); ++i) {
            file << model.pressure(state.A[el.dofOffset() + i], p) << "\n";
        }
    }

    file << "SCALARS TotalPressure double 1\nLOOKUP_TABLE default\n";
    for (const dg::Element& el : elements) {
        const VesselParameters& p = el.vesselParameters();
        for (Index i = 0; i < el.numDofs(); ++i) {
            const Index idx = el.dofOffset() + i;
            file << model.totalPressure({state.A[idx], state.Q[idx]}, p) << "\n";
        }
    }

    file << "SCALARS VesselId int 1\nLOOKUP_TABLE default\n";
    for (const dg::Element& el : elements) {
        for (Index i = 0; i < el.numDofs(); ++i) file << el.vesselId() << "\n";
    }
}

VtkSeriesWriter::VtkSeriesWriter(std::filesystem::path directory, std::string baseName)
    : directory_(std::move(directory)), baseName_(std::move(baseName)) {
    std::filesystem::create_directories(directory_);
}

void VtkSeriesWriter::write(const dg::Mesh& mesh, const physics::State& state,
                             const physics::BloodFlowModel& model, Real time) {
    std::ostringstream nameStream;
    nameStream << baseName_ << '_' << std::setfill('0') << std::setw(5) << nextIndex_ << ".vtk";
    const std::string fileName = nameStream.str();

    writeVtk(directory_ / fileName, mesh, state, model, time);
    writtenSnapshots_.emplace_back(nextIndex_, time);
    ++nextIndex_;

    writeCollectionFile();
}

void VtkSeriesWriter::writeCollectionFile() const {
    std::ofstream file(directory_ / (baseName_ + ".pvd"));
    if (!file) {
        throw std::runtime_error("VtkSeriesWriter: cannot open collection file for '" + baseName_ + "'");
    }
    file << std::setprecision(10);
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    file << "  <Collection>\n";
    for (const auto& [index, time] : writtenSnapshots_) {
        std::ostringstream nameStream;
        nameStream << baseName_ << '_' << std::setfill('0') << std::setw(5) << index << ".vtk";
        file << "    <DataSet timestep=\"" << time << "\" file=\"" << nameStream.str() << "\"/>\n";
    }
    file << "  </Collection>\n";
    file << "</VTKFile>\n";
}

} // namespace hemo1d::output