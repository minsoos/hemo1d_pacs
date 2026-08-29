#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include <catch2/catch_test_macros.hpp>
 
#include "hemo1d/dg/mesh.hpp"
#include "hemo1d/output/vtk_writer.hpp"

using namespace hemo1d;
using namespace hemo1d::dg;
using namespace hemo1d::output;
using namespace hemo1d::physics;

namespace {

Network makeSingleVesselNetwork() {
    VesselParameters params;
    params.length = 2.0;
    params.A0 = 0.126;
    params.beta = 6.06e5;
    params.nElements = 3;
    params.polynomialOrder = 2;

    std::vector<Vessel> vessels{Vessel(1, "v", params)};
    BoundaryConditionSpec bc;
    bc.type = BoundaryConditionType::NonReflecting;
    std::vector<Node> nodes{
        Node(1, "in", {{1, VesselEnd::Proximal}}, bc),
        Node(2, "out", {{1, VesselEnd::Distal}}, bc),
    };
    return Network(FluidProperties{}, std::move(vessels), std::move(nodes));
}
 
std::string readFile(const std::filesystem::path& path) {
    std::ifstream file(path);
    std::ostringstream ss;
    ss << file.rdbuf();
    return ss.str();
}
 
} // namespace


TEST_CASE("writeVtk produces a well-formed legacy VTK PolyData file", "[vtk_writer]") {
    const Network network = makeSingleVesselNetwork();
    const Mesh mesh(network);

    State state(mesh.totalDofs());
    for (Index i = 0; i < state.size(); ++i) {
        state.A[i] = 0.126;
        state.Q[i] = 0.05;
    }

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "hemo1d_test_snapshot.vtk";
    LinearElasticTubeLaw law;
    writeVtk(path, mesh, state, law, 1.05, 0.5);

    REQUIRE(std::filesystem::exists(path));

    const std::string content = readFile(path);
    CHECK(content.find("# vtk DataFile Version 3.0") != std::string::npos);
    CHECK(content.find("DATASET POLYDATA") != std::string::npos);

    const Index expectedPoints = mesh.totalDofs();
    CHECK(content.find("POINTS " + std::to_string(expectedPoints) + " double") != std::string::npos);
    CHECK(content.find("LINES " + std::to_string(mesh.elementCount())) != std::string::npos);
    CHECK(content.find("POINT_DATA " + std::to_string(expectedPoints)) != std::string::npos);
    CHECK(content.find("SCALARS Area double 1") != std::string::npos);
    CHECK(content.find("SCALARS FlowRate double 1") != std::string::npos);
    CHECK(content.find("SCALARS Pressure double 1") != std::string::npos);
    CHECK(content.find("SCALARS TotalPressure double 1") != std::string::npos);
    CHECK(content.find("SCALARS Velocity double 1") != std::string::npos);
    CHECK(content.find("SCALARS VesselId int 1") != std::string::npos);
}

TEST_CASE("VtkSeriesWriter writes numbered snapshots and a pvd collection", "[vtk_writer]") {
    const Network network = makeSingleVesselNetwork();
    const Mesh mesh(network);

    State state(mesh.totalDofs());
    for (Index i = 0; i < state.size(); ++i) {
        state.A[i] = 0.126;
        state.Q[i] = 0.0;
    }

    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "hemo1d_test_series";
    LinearElasticTubeLaw law;
 
    VtkSeriesWriter writer(dir, "series");
    writer.write(mesh, state, law, 1.05, 0.0);

    for (Index i = 0; i < state.size(); ++i) {
        state.A[i] +=0.01;
        state.Q[i] +=0.01;
    }

    writer.write(mesh, state, law, 1.05, 0.1);

    REQUIRE(std::filesystem::exists(dir / "series_00000.vtk"));
    REQUIRE(std::filesystem::exists(dir / "series_00001.vtk"));
    REQUIRE(std::filesystem::exists(dir / "series.pvd"));

    const std::string pvd = readFile(dir / "series.pvd");
    CHECK(pvd.find("series_00000.vtk") != std::string::npos);
    CHECK(pvd.find("series_00001.vtk") != std::string::npos);
    CHECK(pvd.find("<Collection>") != std::string::npos);
}