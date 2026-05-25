# Validation and Testing

OWENSAero validation has two layers: pinned regression tests for code behavior and physics-facing example comparisons. Both are needed because several public workflows are stateful and model options interact.

## Current Test Evidence

| Area | Evidence |
| --- | --- |
| Slice DMS and AC | `test/simple_example.jl` pins output tuples, wind-direction cases, and HDF5 fixture comparisons. |
| Boeing-Vertol dynamic stall | `test/dyn_stall_tests.jl` pins the full NACA 0012 pitching histories and representative extrema without GUI plotting. |
| Unsteady method | `test/simple_example.jl` covers RPI and non-RPI cases, including InflowWind-driven branches. |
| Gradients | steady DMS gradient checks compare AD behavior for representative calls; AC gradient coverage remains open. |
| Full turbine | `test/full_turb.jl` and `test/full_turb_undersampling.jl` exercise the stateful turbine workflow. |
| API sanity | `test/api_unit_tests.jl` pins constructors, types, selected helper behavior, AeroDyn blade/polar readers, the HAWT CCBlade adapter, and the standalone Prandtl tip/root-loss primitive. |

See [Aero Model Audit](@ref) for the current coupled-versus-helper status of each aerodynamic addition.

## Validation Rules for New Cases

- Pin concrete numeric baselines with absolute or relative tolerances tied to model accuracy, not just `isreal` or nonzero checks.
- Record model choices in the test name: `AeroModel`, `DynamicStallModel`, `RPI`, `tau`, InflowWind state, added-mass state, and buoyancy state.
- Check output shape, element type, and at least representative force, velocity, angle, and coefficient values.
- For coupled OWENS tests, assert transformed loads after frame mapping as well as local aerodynamic loads.
- Add gradient checks for any branch intended for optimization.

## Open Validation Gaps

The RM2 and SNL example scripts are useful physics references, but not every comparison is currently a CI-quality validation. Remaining hardening work includes expanding added-mass and buoyancy baselines beyond aggregate values, normalizing the HAWT OpenFAST/AeroDyn root-station and torque/power frame conventions before comparing against CCBlade, and pinning mode/coupled-load transformations in OWENS integration tests.
