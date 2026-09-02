#include <functional>
#include <utility>
#include <vector>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "hemo1d/io/network_parser.hpp"
#include "hemo1d/couplings/register.hpp"
#include "hemo1d/couplings/windkessel_coupling.hpp"
#include "hemo1d/physics/section_state.hpp"
#include "hemo1d/physics/terminal_coupling.hpp"
#include "hemo1d/simulation.hpp"

namespace py = pybind11;
using namespace hemo1d;

PYBIND11_MODULE(hemo1d, m) {
    m.doc() = "Hemo1D: a 1D hemodynamics Discontinuous Galerkin solver";

    couplings::registerBuiltinCouplings();

    py::enum_<FluxKind>(m, "FluxKind")
        .value("LAX_FRIEDRICHS", FluxKind::LaxFriedrichs)
        .value("HLL", FluxKind::Hll);

    py::class_<Network>(m, "Network")
        .def_property_readonly("vessel_count", &Network::vesselCount)
        .def_property_readonly("node_count", &Network::nodeCount);

    m.def(
        "load_network", [](const std::string& path) { return io::loadNetwork(path); }, py::arg("path"),
        "Load a vessel network topology from a JSON file (see examples/networks/).");

    py::class_<SimulationSettings>(m, "SimulationSettings")
        .def(py::init<>())
        .def_readwrite("default_polynomial_order", &SimulationSettings::defaultPolynomialOrder)
        .def_readwrite("flux", &SimulationSettings::flux)
        .def_readwrite("use_slope_limiter", &SimulationSettings::useSlopeLimiter)
        .def_readwrite("slope_limiter_tolerance", &SimulationSettings::slopeLimiterTolerance);

    py::class_<physics::State>(m, "State")
        .def_readonly("A", &physics::State::A)
        .def_readonly("Q", &physics::State::Q);

    // --- exterior coupling (0D / external systems at terminal nodes) ---------

    py::class_<physics::TerminalInterface>(m, "TerminalInterface", "Interior interface state handed to a coupling callback.")
        .def_property_readonly(
            "end", [](const physics::TerminalInterface& i) { return i.end == VesselEnd::Proximal ? "proximal" : "distal"; })
        .def_property_readonly("A", [](const physics::TerminalInterface& i) { return i.trace.A; })
        .def_property_readonly("Q", [](const physics::TerminalInterface& i) { return i.trace.Q; })
        .def_property_readonly("dAdz", [](const physics::TerminalInterface& i) { return i.grad.dAdz; })
        .def_property_readonly("dQdz", [](const physics::TerminalInterface& i) { return i.grad.dQdz; })
        .def_property_readonly("rho", [](const physics::TerminalInterface& i) { return i.rho; })
        .def_property_readonly("A0", [](const physics::TerminalInterface& i) { return i.params.A0; })
        .def_property_readonly("beta", [](const physics::TerminalInterface& i) { return i.params.beta; })
        .def_property_readonly("alpha", [](const physics::TerminalInterface& i) { return i.params.alpha; })
        .def_property_readonly("length",
                                [](const physics::TerminalInterface& i) { return i.params.length; })
        .def_property_readonly("friction_kr",
                                [](const physics::TerminalInterface& i) { return i.params.frictionKr; });

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
        .def(py::init<Network, SimulationSettings>(), py::arg("network"),
             py::arg("settings") = SimulationSettings{})

        .def("set_uniform_initial_condition", &Simulation::setUniformInitialCondition,
             py::arg("flow_rate") = 0.0, py::arg("area") = -1.0,
             "Sets every DOF to the given flow rate and to its vessel's A0 "
             "(or to `area` uniformly, if positive).")

        .def(
            "set_windkessel_outlet",
            [](Simulation& sim, Id nodeId, Real r1,
               std::vector<std::pair<Real, Real>> compartments, Real pOut, Real pInit, int subSteps,
               const std::string& pOutCsv) {
                couplings::WindkesselParameters wk;
                wk.R1 = r1;
                for (const auto& rc : compartments) {
                    wk.compartments.push_back(couplings::LumpedCompartment{rc.first, rc.second});
                }
                wk.pOut = pOut;
                wk.pOutCsv = pOutCsv;
                wk.pInit = pInit;
                wk.subSteps = subSteps;
                sim.setCoupling(nodeId, couplings::makeWindkesselCoupling(std::move(wk)));
            },
            py::arg("node_id"), py::arg("r1") = -1.0,
            py::arg("compartments") = std::vector<std::pair<Real, Real>>{}, py::arg("p_out") = 0.0,
            py::arg("p_init") = 0.0, py::arg("sub_steps") = 1, py::arg("p_out_csv") = std::string(),
            "Attach an RC-ladder Windkessel 0D model to a terminal node, replacing "
            "its network boundary condition. `compartments` is a list of (R, C) "
            "pairs (>= 1); r1 < 0 uses the matched impedance rho*c0/A0.")

        .def(
            "set_coupling_callback",
            [](Simulation& sim, Id nodeId, py::function fn, py::object commitFn) {
                Simulation::CouplingCommitFn commit;
                if (!commitFn.is_none()) {
                    py::function cfn = commitFn.cast<py::function>();
                    commit = [cfn](physics::SectionState resolved, const physics::TerminalInterface& iface,
                                   Real t, Real dt) {
                        py::gil_scoped_acquire gil;
                        cfn(std::make_pair(resolved.A, resolved.Q), iface, t, dt);
                    };
                }
                sim.setCouplingCallback(
                    nodeId,
                    [fn](const physics::TerminalInterface& iface, Real t,
                         Real dt) -> physics::SectionState {
                        py::gil_scoped_acquire gil;
                        const auto rc = fn(iface, t, dt).cast<std::pair<Real, Real>>();
                        return {rc.first, rc.second};
                    },
                    std::move(commit));
            },
            py::arg("node_id"), py::arg("fn"), py::arg("commit_fn") = py::none(),
            "Attach a per-step feedback law to a terminal node -- the generic hook "
            "for coupling a 3D subdomain or any external model. "
            "fn(iface, time, dt) -> (A, Q) returns the ghost state; the optional "
            "commit_fn((A, Q), iface, time, dt) is called after the accepted step "
            "so an external model can advance its own state. `iface` is a "
            "TerminalInterface.")

        .def("coupling_state", &Simulation::couplingState, py::arg("node_id"),
             "Scalar internal state of a terminal coupling (e.g. Windkessel "
             "compartment pressures, interface-first); empty for a stateless BC.")

        .def("add_probe", &Simulation::addProbe, py::arg("name"), py::arg("vessel_id"), py::arg("z"),
             py::return_value_policy::reference_internal)

        .def("probe_samples", &Simulation::probeSamples, py::arg("name"))

        .def("write_probes_csv", &Simulation::writeProbesCsv, py::arg("directory"))

        .def("field_snapshot", &Simulation::fieldSnapshot,
             "Returns every degree of freedom's (element_index, vessel_id, z, area, flow_rate) "
             "at the current state.")

        .def("step", &Simulation::step, py::arg("dt"), "Advances the simulation by one step of size dt.")

        .def("cfl_time_step", &Simulation::cflTimeStep, py::arg("cfl_number") = 0.9,
             "Largest stable explicit dt for the current state (DG CFL bound "
             "h/((2p+1)*wave_speed), minimized over elements), scaled by "
             "cfl_number. Recompute after the state changes materially.")

        .def("run", &Simulation::run, py::arg("target_time"), py::arg("dt"), py::arg("record_every") = 1,
             py::arg("vtk_directory") = std::string(), py::arg("vtk_every") = 0,
             py::arg("cfl_number") = 0.9,
             "Advances until target_time, recording probes every record_every "
             "steps and (if vtk_directory is given) writing a VTK snapshot "
             "every vtk_every steps. If dt > 0 it is used fixed for every "
             "step; if dt <= 0 a dynamic CFL-based dt is recomputed every "
             "step instead (see cfl_time_step()), using cfl_number.")

        .def_property_readonly("time", &Simulation::time)

        .def_property_readonly("state", &Simulation::state, py::return_value_policy::reference_internal)

        .def_property_readonly("network", &Simulation::network, py::return_value_policy::reference_internal);
}
