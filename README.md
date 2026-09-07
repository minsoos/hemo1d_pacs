# Hemo1D

Hemo1D is a C++20 solver for one-dimensional blood-flow simulations in vascular networks. The spatial discretization is based on the discontinuous Galerkin method, with Python bindings for setting up simulations, running them, and collecting results.

The repository currently provides:

- vascular networks defined from JSON files;
- HLL and Lax-Friedrichs numerical fluxes;
- prescribed and non-reflecting boundary conditions;
- vessel junctions;
- Windkessel terminal models;
- user-defined terminal coupling from Python;
- fixed or CFL-controlled time stepping;
- probe output to CSV and field output to VTK.

## Installation

Python 3.8 or newer and a C++20-capable compiler are required.

From the repository root:

```bash
python -m pip install .
```

For development, it is usually convenient to use a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
```

Some examples and analysis scripts additionally use NumPy, Pandas, and Matplotlib. Install them as needed:

```bash
python -m pip install numpy pandas matplotlib
```

The first build may download C++ dependencies through CMake.

## Quick start

The simplest way to run the solver is through the Python interface.

```python
import hemo1d

network = hemo1d.load_network(
    "examples/networks/simple_bifurcation_sine.json"
)

settings = hemo1d.SimulationSettings()
settings.default_polynomial_order = 1
settings.flux = hemo1d.FluxKind.HLL
settings.use_slope_limiter = True

sim = hemo1d.Simulation(network, settings)

sim.add_probe("inlet", vessel_id=1, z=0.05)

sim.run(
    target_time=0.01,
    dt=-1.0,              # use a CFL-based time step
    cfl_number=0.8,
    record_every=10,
)

sim.write_probes_csv("out/probes")
```

See [`examples/`](examples/) for complete runnable examples.

## Network input

Networks are described by JSON files containing:

- fluid properties;
- vessels and their discretization;
- network nodes and vessel connectivity;
- inlet and outlet boundary conditions.

Prescribed time-dependent data are read from CSV files.

Example networks are available in [`examples/networks/`](examples/networks/). See its README for the input format and unit conventions.

## Running the examples

For example:

```bash
python examples/run_simple_bifurcation.py
python examples/run_windkessel_rcr.py
python examples/run_python_windkessel.py
```

The examples demonstrate prescribed inlet data, bifurcations, Windkessel outlets, probes, VTK output, and Python-defined coupling.

## Output

Probe values can be recorded at arbitrary positions along vessels and written to CSV:

```python
sim.add_probe("probe_name", vessel_id=1, z=0.5)
sim.write_probes_csv("out/probes")
```

Field snapshots can also be written as VTK files by passing a directory to `Simulation.run`:

```python
sim.run(
    target_time=0.1,
    dt=-1.0,
    vtk_directory="out/vtk",
    vtk_every=100,
)
```

The generated VTK series can be opened with ParaView.

## Units

The examples use a consistent CGS unit system:

| Quantity | Unit |
|---|---|
| length | cm |
| area | cm² |
| time | s |
| velocity | cm/s |
| flow rate | cm³/s |
| density | g/cm³ |
| dynamic viscosity | g/(cm s) |
| pressure | dyn/cm² |

Network and boundary-condition data must use mutually consistent units.

## C++ build

To build the C++ library and tests directly:

```bash
cmake -S . -B build -DHEMO1D_BUILD_PYTHON=OFF
cmake --build build -j
ctest --test-dir build --output-on-failure
```

CMake 3.20 or newer is required.

A small benchmark executable is also built:

```bash
./build/apps/hemo1d_bench
```

An optional integer argument controls the number of benchmark repetitions.

## Repository layout

| Path | Purpose |
|---|---|
| `include/hemo1d/` | public C++ headers |
| `src/` | C++ implementation |
| `python/` | Python bindings |
| `examples/` | small runnable examples |
| `examples/networks/` | example network and inlet files |
| `tests/` | C++ unit tests |
| `convergence/` | numerical convergence studies |
| `experiments/CoW/` | Circle of Willis experiment |
| `apps/` | utility and benchmark programs |

For numerical verification, see [`convergence/`](convergence/). For the larger Circle of Willis case, see [`experiments/CoW/`](experiments/CoW/).
