# Model Selection

OWENSAero exposes several modeling choices. Select them explicitly because each option changes both the physics represented and the test evidence that applies to a run.

## DMS

Set `AeroModel = "DMS"` for the double-multiple-streamtube path. This is the default high-throughput model used in many OWENS examples and optimization-style tests. It solves azimuthal induction and returns sectional force, velocity, angle-of-attack, coefficient, and Reynolds-number arrays.

Use DMS when:

- the run needs fast repeated evaluations;
- automatic differentiation matters;
- the turbine is represented as stacked independent vertical slices;
- the requested fidelity is consistent with streamtube assumptions.

## Actuator Cylinder

Set `AeroModel = "AC"` for the actuator-cylinder path. The implementation is based on the BYU FLOW Lab VAWT AC method and has been adapted for OWENS slice and turbine workflows.

Use AC when:

- actuator-cylinder loading is the intended validation basis;
- multi-turbine or radial-influence behavior is being investigated;
- the extra cost is acceptable relative to DMS.

## Finite-Span and Tip-Loss Scope

DMS and AC currently operate as stacked two-dimensional slice models. The
current regression suite pins the no-tip-loss baseline: moving an otherwise
identical slice from mid-span to a tip-adjacent `z` location does not change the
computed CP, torque, or distributed loads. Any future finite-span or tip-loss
correction should introduce an explicit model option and validation data, then
update that baseline with tests that show the expected numerical change.

## Dynamic Stall

`DynamicStallModel = "BV"` activates the Boeing-Vertol model used by the current tests and examples. `DynamicStallModel = "none"` bypasses dynamic stall and is the safest choice for gradient-sensitive workflows. The Leishman-Beddoes path is present in historical comments but should not be treated as a production option until it has complete tests and documentation.

Model option strings are normalized at construction and setup time: `BV`/`boeing-vertol` select Boeing-Vertol, `none`/`noDS`/`NONE` disable dynamic stall, and `DMS`/`AC` select the aerodynamic model case-insensitively. `DynamicStallModel = "LB"` now fails immediately with a not-implemented error instead of falling through to a static airfoil path.

## Unsteady Method and RPI

`UnsteadyParams` controls the unsteady wake path. `RPI = true` enables reduced-pass interpolation to avoid recalculating every azimuthal point at every unsteady step. The `tau` vector controls wake propagation filtering. Existing tests pin representative unsteady results for both DMS and AC, but new use cases should add explicit baseline values for the selected `tau`, `RPI`, and time-step choices.

## InflowWind

Set `ifw = true` in `UnsteadyParams` or the high-level setup when turbulent inflow is provided by OpenFAST InflowWind. The `.bts` files are external TurbSim/OpenFAST inputs; the test data under `test/data/ifw` shows the required layout.

## Added Mass, Buoyancy, and Rotational Acceleration Terms

The high-level `setupTurb` path accepts flags for added-mass, buoyancy, rotational acceleration, and centrifugal terms. These are important for water-turbine or fast-transient studies. They also change the returned load components, so they should be pinned in validation tests before being used as acceptance evidence for a new design.
