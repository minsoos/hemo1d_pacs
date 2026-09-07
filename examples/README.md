# Examples

This directory contains small Python examples illustrating the main Hemo1D workflows.

Install Hemo1D from the repository root before running them:

```bash
python -m pip install .
```

Matplotlib is used for plotting:

```bash
python -m pip install matplotlib
```

## Simple bifurcation

```bash
python examples/run_simple_bifurcation.py
```

Runs a three-vessel bifurcation with a time-dependent prescribed inlet.

The example demonstrates:

- loading a network from JSON;
- configuring the DG solver;
- adding probes;
- running with a fixed time step;
- writing probe CSV files;
- writing VTK snapshots;
- plotting the resulting waveforms.

The default network is:

```text
examples/networks/simple_bifurcation_sine.json
```

## Windkessel outlet

```bash
python examples/run_windkessel_rcr.py
```

Runs a single vessel terminated by an RCR-type Windkessel model.

Besides the vessel solution, the script records the internal Windkessel state during the simulation.

## Python coupling callback

```bash
python examples/run_python_windkessel.py
```

Demonstrates how an external terminal model can be implemented from Python.

The script compares the Python callback implementation with the built-in C++ Windkessel coupling.

## Network files

Input networks and prescribed time series are under:

```text
examples/networks/
```

See [`networks/README.md`](networks/README.md) for the JSON structure, boundary conditions, CSV format, and units.

## Output

Each script writes into its own output directory unless another directory is supplied.

Typical outputs include:

```text
probes/
    <probe>.csv
vtk/
    ...
```

Some examples also produce plots or additional CSV files containing coupling states.
