# Validation and Testing

OWENSAero validation has two layers: pinned regression tests for code behavior and physics-facing example comparisons. Both are needed because several public workflows are stateful and model options interact.

## Current Test Evidence

| Area | Evidence |
| --- | --- |
| Slice DMS and AC | `test/simple_example.jl` pins output tuples, wind-direction cases, and HDF5 fixture comparisons. |
| Boeing-Vertol dynamic stall | `test/dyn_stall_tests.jl` pins the full NACA 0012 pitching histories and representative extrema without GUI plotting. |
| Leishman-Beddoes dynamic stall | `test/dyn_stall_tests.jl` pins NACA 0012 pitching histories, leading-edge separation state, vortex-age state, and branch metrics against digitized Sadr Leishman-Beddoes curves. |
| Unsteady method | `test/simple_example.jl` covers RPI and non-RPI cases, including InflowWind-driven branches. |
| Gradients | steady DMS gradient checks compare AD behavior for representative calls; AC gradient coverage remains open. |
| Full turbine | `test/full_turb.jl` and `test/full_turb_undersampling.jl` exercise the stateful turbine workflow. |
| API sanity | `test/api_unit_tests.jl` pins constructors, types, selected helper behavior, AeroDyn blade/polar readers, the HAWT CCBlade adapter, and the standalone Prandtl tip/root-loss primitive. |
| HAWT AeroDyn BEM | `test/api_unit_tests.jl` runs the native AeroDyn Basic HAWT driver on the NREL 5 MW aerodynamic fixture and compares CCBlade rotor totals plus interior station `Alpha`, `Cl`, `Cd`, axial induction, coefficient, and mapped load channels. |

See [Aero Model Audit](@ref) for the current coupled-versus-helper status of each aerodynamic addition.

## Validation Rules for New Cases

- Pin concrete numeric baselines with absolute or relative tolerances tied to model accuracy, not just `isreal` or nonzero checks.
- Record model choices in the test name: `AeroModel`, `DynamicStallModel`, `RPI`, `tau`, InflowWind state, added-mass state, and buoyancy state.
- Check output shape, element type, and at least representative force, velocity, angle, and coefficient values.
- For coupled OWENS tests, assert transformed loads after frame mapping as well as local aerodynamic loads.
- Add gradient checks for any branch intended for optimization.

## Open Validation Gaps

The RM2 and SNL example scripts are useful physics references, but not every comparison is currently a CI-quality validation. Remaining hardening work includes expanding added-mass and buoyancy baselines beyond aggregate values, extending the HAWT AeroDyn comparison beyond steady Basic HAWT BEM into DBEMT and moving-frame cases, and pinning mode/coupled-load transformations in OWENS integration tests.
