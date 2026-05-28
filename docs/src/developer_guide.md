# Developer Guide

This guide is for changes to OWENSAero models, examples, tests, and docs. Most
public workflows are stateful and feed coupled structural simulations, so small
changes should be reviewed for numerical and physical consistency.

## Local Workflow

Use the tests and examples together:

| Area | Primary files |
|:---|:---|
| Direct DMS/AC slices | `test/simple_example.jl`, `test/dms_unit_tests.jl` |
| Full turbine lifecycle | `test/full_turb.jl`, `test/full_turb_undersampling.jl` |
| Dynamic stall | `test/dyn_stall_tests.jl`, `examples/dynamic_stall` |
| HAWT adapter | `test/hawt_example_tests.jl`, `examples/HAWT/rigid_ccblade_hawt.jl` |
| AeroDyn readers | `test/api_unit_tests.jl`, `src/AeroDynInputs.jl` |

Run focused tests while developing, then run the package test suite before
opening a PR:

```sh
julia --project=. test/runtests.jl
```

Build docs from the package root:

```sh
julia --project=docs docs/make.jl
```

## Model Changes

Before changing DMS, AC, or full-turbine code, check:

| Topic | Review point |
|:---|:---|
| Units | Keep all public inputs and outputs in SI units. Angles are radians unless a name explicitly carries `_d`. |
| Frames | Document any change to radial, tangential, vertical, blade, hub, or global load directions. |
| Signs | Check torque, `CP`, tangential force, wind-angle rotation, clockwise/counter-clockwise cases, and pitching moments independently. |
| State | Reinitialize with `setupTurb` when topology, airfoil tables, inflow source, or discretization changes. |
| Tuple contracts | Update tests and docs when steady or unsteady output tuple order changes. |
| Interpolation | Keep airfoil interpolation inside validated tabular bounds unless the extrapolation behavior is explicitly tested. |
| AD | Add derivative checks for branches expected to support optimization. |

## Test Expectations

Tests should pin concrete values, not just finite or nonzero results. For a new
branch, include at least one output-shape assertion, one force or moment value,
one local velocity or angle-of-attack value, and one coefficient or state value.

For coupled OWENS-facing changes, test the transformed load channel as well as
the local aerodynamic value. This catches sign and ordering mistakes that are
invisible in aggregate power.

When adding examples, keep generated output out of git unless the file is a
small reviewed fixture needed by a deterministic test.

## Documentation

Docs should describe the model that exists in the current checkout. If a page
mentions a validation case, point to the test or example that exercises it.
Public docstrings should name mutation, state ownership, units, and output
ordering whenever those details matter for coupled simulations.
