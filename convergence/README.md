# Convergence studies

This directory contains scripts used to study numerical convergence of Hemo1D.

These scripts are intended for solver verification rather than as general-purpose simulation examples.

## Requirements

Install Hemo1D and the analysis dependencies:

```bash
python -m pip install .
python -m pip install numpy matplotlib
```

## Files

`common.py`
: Shared physical parameters, mesh sizes, polynomial orders, and helper functions.

`run_convergence.py`
: Generates network inputs and runs the convergence cases.

`analyze_convergence.py`
: Compares the generated solutions across resolutions and estimates convergence rates.

Generated files are written below:

```text
convergence/generated/
convergence/output/
```

## Running the spatial study

From the repository root:

```bash
python convergence/run_convergence.py
```

The current `main()` runs the single-vessel spatial convergence study.

For each resolution and polynomial order, the script stores simulation metadata and field snapshots.

Typical case output is:

```text
manifest.json
field_snapshots.csv
```

## Other studies

`run_convergence.py` also contains functions for bifurcation and Windkessel convergence cases:

```python
run_bifurcation_cases(...)
run_windkessel_cases(...)
```

Their calls are currently disabled in `main()`. Enable the corresponding call when those data are needed.

## Analysis

Analyze the spatial study with:

```bash
python convergence/analyze_convergence.py \
    --study spatial \
    --quantity flow_rate
```

Available studies are:

```text
spatial
bifurcation
windkessel
```

The quantity can be:

```text
area
flow_rate
```

For example:

```bash
python convergence/analyze_convergence.py \
    --study spatial \
    --quantity area
```

The analysis script uses Richardson extrapolation by default and writes the resulting convergence plot below `convergence/output/`.

Use:

```bash
python convergence/analyze_convergence.py --help
```

for the complete set of options.

## Changing the study

Mesh sizes, polynomial orders, time steps, physical parameters, and target times are defined in `common.py`.

When changing these values, keep the time step sufficiently small that temporal error does not dominate a spatial-convergence study.
