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
- single-turbine radial-influence behavior is being investigated;
- the extra cost is acceptable relative to DMS.

The internal actuator-cylinder implementation contains multi-turbine assembly
building blocks, but the public `AC` path currently rejects multi-turbine calls.
Treat multi-turbine AC as open solver work until the public path is restored
and validated.

## Finite-Span and Tip-Loss Scope

DMS and AC currently operate as stacked two-dimensional slice models. The
current regression suite pins the no-tip-loss baseline: moving an otherwise
identical slice from mid-span to a tip-adjacent `z` location does not change the
computed CP, torque, or distributed loads. Any future finite-span or tip-loss
correction should introduce an explicit model option and validation data, then
update that baseline with tests that show the expected numerical change.

`prandtlTipLossFactor(num_blades, radial_position, rotor_radius, inflow_angle;
hub_radius=0.0, include_root=false)` provides a pinned Prandtl finite-blade
loss primitive for caller-side axial-flow studies and future HAWT work. It is
not automatically applied inside DMS or AC because those VAWT solvers need a
separate finite-height validation basis before their load baselines should
change.

## HAWT Dynamic Inflow Scope

`oyeDynamicInflowTimeConstants`, `oyeDynamicInflowDerivative`, and
`oyeDynamicInflowStep` provide a tested Oye dynamic-inflow primitive for future
CCBlade-based HAWT work. The helper follows the AeroDyn DBEMT continuous
state-space form and advances axial or tangential induction states at radial
stations. It is not a replacement for the current DMS wake-speed filter and is
not yet coupled into a HAWT solver path; callers must provide the quasi-steady
induction from a BEM backend such as CCBlade and map the resulting dynamic
induction into loads explicitly.

## Dynamic Stall

`DynamicStallModel = "BV"` activates the Boeing-Vertol model used by the current tests and examples. `DynamicStallModel = "none"` bypasses dynamic stall and is the safest choice for gradient-sensitive workflows. The Leishman-Beddoes path is present in historical comments but should not be treated as a production option until it has complete tests and documentation.

Model option strings are normalized at construction and setup time: `BV`/`boeing-vertol` select Boeing-Vertol, `none`/`noDS`/`NONE` disable dynamic stall, and `DMS`/`AC` select the aerodynamic model case-insensitively. `DynamicStallModel = "LB"` now fails immediately with a not-implemented error instead of falling through to a static airfoil path.

## Parasitic and Lumped Drag Scope

`jointDragForce(rho, velocity, CdA)` provides a pinned lumped bluff-body drag primitive for joints, hardware, or other ancillary bodies. It returns a force vector in the same frame as the supplied relative velocity using `-0.5*rho*CdA*|V|*V`. This helper is intentionally not coupled into DMS or AC induction, and it does not modify blade airfoil loads or turbine base loads unless a caller explicitly maps the returned force into a structural load workflow.

## Tower Shadow Scope

`towerShadowVelocity(velocity, relative_position; active=false, ...)` provides a pinned, off-by-default inflow-deficit primitive for tower shadow studies. The active model applies a Gaussian downstream velocity deficit in the supplied wind frame. It is intentionally not coupled into DMS or AC induction yet; callers must explicitly use the returned velocity in their own load-building workflow until a validated solver-integrated tower-shadow option is added.

## Lifting Strut Scope

`liftingStrutForce(rho, velocity, chord, span, cl, cd, lift_direction)` provides a pinned integrated segment-force primitive for lifting struts or other non-blade lifting members. It returns `0.5*rho*|V|^2*chord*span*(cl*lift_hat - cd*velocity_hat)`, where drag opposes the supplied relative velocity and the caller supplies the strut-geometry-projected positive-lift direction in the same frame. This helper is intentionally not coupled into DMS or AC induction yet; callers must explicitly map the returned force into the structural load workflow until a validated solver-integrated strut option is added.

## Unsteady Method and RPI

`UnsteadyParams` controls the unsteady wake path. `RPI = true` enables reduced-pass interpolation to avoid recalculating every azimuthal point at every unsteady step. The `tau` vector controls wake propagation filtering. Existing tests pin representative unsteady results for both DMS and AC, but new use cases should add explicit baseline values for the selected `tau`, `RPI`, and time-step choices.

`tau` must contain two finite positive values. The filtered wake update uses
the magnitude of the prior wake speed with a small positive floor, so a
transient negative or zero wake estimate does not create zero time constants.
RPI with negative mean rotor speed is intentionally rejected because the
clockwise/negative-RPM unsteady-RPI path is not yet validated; use `RPI=false`
for negative-RPM regression work until the full clockwise unsteady path is
implemented and pinned.

## InflowWind

Set `ifw = true` in `UnsteadyParams` or the high-level setup when turbulent inflow is provided by OpenFAST InflowWind. The `.bts` files are external TurbSim/OpenFAST inputs; the test data under `test/data/ifw` shows the required layout.

## Added Mass, Buoyancy, and Rotational Acceleration Terms

The high-level `setupTurb` path accepts flags for added-mass, buoyancy, rotational acceleration, and centrifugal terms. These are important for water-turbine or fast-transient studies. They also change the returned load components, so they should be pinned in validation tests before being used as acceptance evidence for a new design.
