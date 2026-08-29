#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hemo1d/io/network_parser.hpp"
#include "hemo1d/simulation.hpp"

namespace py = pybind11;
using namespace hemo1d;

PYBIND11_MODULE(hemo1d, m) {
    m.doc() = "Hemo1d: a 1D hemodynamics Dicontinuous Galerkin solver";

    py::enum_<FluxKind>(m, "FluxKind")
        .value("LAX_FRIEDRICHS", FluxKind::LaxFriedrichs)
        .value("HLL", FluxKind::Hll);

    py::class_<Network>(m, "Network")
        .def_property_readonly("vessel_count", &Network::vesselCount)
        .def_property_readonly("node_count", &Network::nodeCount);
    
    
    m.def(
        "load_network", [](const std::string& path) { return io::loadNetwork(path); }, py::arg("path"),
        "Load a vessel network topology from a JSON file."
    );

    py::class_<SimulationSettings>(m, "SimulationSettings")
        .def(py::init<>())
        .def_readwrite("default_polynomial_order", &SimulationSettings::defaultPolynomialOrder)
        .def_readwrite("flux", &SimulationSettings::flux)
        .def_readwrite("use_slope_limiter", &SimulationSettings::useSlopeLimiter)
        .def_readwrite("slope_limiter_tolerance", &SimulationSettings::slopeLimiterTolerance);

    py::class_<physics::State>(m, "State")
        .def_readonly("A", &physics::State::A)
        .def_readonly("Q", &physics::State::Q);

    py::class_<output::ProbeSample>(m, "ProbeSample")
        .def_readonly("time", &output::ProbeSample::time)
        .def_readonly("area", &output::ProbeSample::A)
        .def_readonly("flow_rate", &output::ProbeSample::Q)
        .def_readonly("pressure", &output::ProbeSample::pressure)
        .def_readonly("velocity", &output::ProbeSample::velocity);

    py::class_<output::Probe>(m, "Probe")
        .def_property_readonly("name", &output::Probe::name)
        .def_property_readonly("vessel_id", &output::Probe::vesselId)
        .def_property_readonly("z", &output::Probe::z);

    py::class_<output::FieldSample>(m, "FieldSample")
        .def_readonly("element_index", &output::FieldSample::elementIndex)
        .def_readonly("vessel_id", &output::FieldSample::vesselId)
        .def_readonly("z", &output::FieldSample::z)
        .def_readonly("area", &output::FieldSample::A)
        .def_readonly("flow_rate", &output::FieldSample::Q);

    py::class_<Simulation>(m, "Simulation")
        .def(py::init<Network, SimulationSettings>(), py::arg("network"), py::arg("settings") = SimulationSettings{})
        
        .def(
            "set_uniform_initial_condition", &Simulation::setUniformInitialCondition,
            py::arg("flow_rate") = 0.0, py::arg("area") = -1.0,
            "Sets every DOF to the given flow rate and to its vessel's A0 (or to `area` uniformly, if positive)."
        )

        .def(
            "add_probe", &Simulation::addProbe,
            py::arg("name"), py::arg("vessel_id"), py::arg("z"),
            py::return_value_policy::reference_internal
        )

        .def(
            "probe_samples", &Simulation::probeSamples,
            py::arg("name")
        )

        .def(
            "write_probes_csv", &Simulation::writeProbesCsv,
            py::arg("directory")
        )

        .def(
            "field_snapshot", &Simulation::fieldSnapshot,
            "Returns every DOF's state value (element_index, vessel_id, z, area, flow_rate) "
            "at the current state"
        )

        .def(
            "step", &Simulation::step,
            py::arg("dt"),
            "Advances the simulation by one step of size dt."
        )

        .def(
            "cfl_time_step", &Simulation::cflTimeStep,
            py::arg("cfl_number") = 0.9,
            "Largest stable value of dt for the current state."
        )

        .def(
            "run", &Simulation::run,
            py::arg("target_time"), py::arg("dt"), py::arg("record_every") = 1,
            py::arg("vtk_directory") = std::string(), py::arg("vtk_every") = 0,
            py::arg("cfl_number") = 0.9,
            "Advances the simulation until `target_time`, recording probes every `record_every` steps "
            "and (if `vtk_directory` is giving) writing a VTK snapshot every "
            "`vtk_every` steps. If dt <= 0, a dynamic CFL-based dt is used."
        )

        .def_property_readonly("time", &Simulation::time)
        .def_property_readonly("state", &Simulation::state, py::return_value_policy::reference_internal)
        .def_property_readonly("network", &Simulation::network, py::return_value_policy::reference_internal);

}
