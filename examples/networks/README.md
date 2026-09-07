# Network input files

Hemo1D vascular networks are defined using JSON files.

A network contains fluid properties, vessels, nodes, and boundary conditions. Time-dependent prescribed data are stored separately as CSV files.

## Minimal structure

A network has the following general form:

```json
{
  "fluid": {
    "density": 1.055,
    "viscosity": 0.045
  },
  "vessels": [
    {
      "id": 1,
      "name": "vessel",
      "length": 10.0,
      "A0": 1.0,
      "beta": 1000000.0,
      "n_elements": 20
    }
  ],
  "nodes": [
    {
      "id": 1,
      "connections": [
        {"vessel": 1, "end": "proximal"}
      ],
      "boundary_condition": {
        "type": "prescribed",
        "quantity": "flow_rate",
        "csv_file": "data/inlet.csv"
      }
    },
    {
      "id": 2,
      "connections": [
        {"vessel": 1, "end": "distal"}
      ],
      "boundary_condition": {
        "type": "non_reflecting"
      }
    }
  ]
}
```

The optional top-level `_description` field can be used to describe the case.

## Fluid

```json
"fluid": {
  "density": 1.055,
  "viscosity": 0.045
}
```

The examples use density in `g/cm³` and dynamic viscosity in `g/(cm s)`.

## Vessels

Each vessel requires:

| Field | Description |
|---|---|
| `id` | unique integer vessel identifier |
| `length` | vessel length |
| `A0` | reference cross-sectional area |
| `beta` | wall stiffness parameter |

Optional fields are:

| Field | Description |
|---|---|
| `name` | human-readable vessel name |
| `alpha` | momentum correction factor; default is `4/3` |
| `friction_kr` | friction coefficient; default is `0` |
| `n_elements` | number of DG elements; default is `1` |
| `polynomial_order` | polynomial order for this vessel |

`length`, `A0`, and `beta` must be positive. `n_elements` and `polynomial_order` must be at least one.

## Nodes and connections

A node lists the vessel ends that meet at that location.

```json
"connections": [
  {"vessel": 1, "end": "distal"},
  {"vessel": 2, "end": "proximal"},
  {"vessel": 3, "end": "proximal"}
]
```

`end` must be either:

```text
proximal
distal
```

Terminal nodes normally have one connection. Junction nodes contain multiple connections.

Bifurcations may also provide:

```json
"bifurcation_angles_rad": [0.7, 1.0]
```

with angles expressed in radians.

## Boundary conditions

### Prescribed value

A prescribed boundary condition reads a time series from CSV:

```json
"boundary_condition": {
  "type": "prescribed",
  "quantity": "flow_rate",
  "csv_file": "data/inlet.csv"
}
```

Supported quantities are:

```text
flow_rate
velocity
area
pressure
```

CSV paths are resolved relative to the network JSON file.

### Non-reflecting outlet

```json
"boundary_condition": {
  "type": "non_reflecting"
}
```

This is useful for simple cases where reflected waves from the terminal boundary should be minimized.

### Windkessel outlet

Windkessel models are specified as external boundary conditions:

```json
"boundary_condition": {
  "type": "external",
  "model": "windkessel",
  "params": {
    "r1": -1.0,
    "compartments": [
      {
        "r": 50000.0,
        "c": 8e-7
      }
    ],
    "p_out": 0.0,
    "p_init": 0.0,
    "sub_steps": 1
  }
}
```

A negative `r1` requests the matched characteristic impedance at the vessel outlet.

See `single_vessel_windkessel.json` and `windkessel_rcr.json` for complete examples.

## Prescribed CSV files

Time-dependent input files contain two columns:

```text
time,value
0.000,0.0
0.010,5.0
0.020,8.0
```

A header row is optional.

Blank lines and lines beginning with `#` are ignored. Sample times must be strictly increasing.

Values between samples are linearly interpolated. Before the first sample and after the last sample, the corresponding endpoint value is used.

## Units

The supplied examples use CGS units:

| Quantity | Unit |
|---|---|
| vessel length | cm |
| area | cm² |
| time | s |
| velocity | cm/s |
| flow rate | cm³/s |
| density | g/cm³ |
| viscosity | g/(cm s) |
| pressure | dyn/cm² |

All parameters in a network should use a consistent unit system.

## Validation

The JSON parser is intentionally strict. Unknown fields and invalid values produce an error rather than being silently ignored.

When creating a new case, starting from one of the files in this directory is usually the easiest approach.
