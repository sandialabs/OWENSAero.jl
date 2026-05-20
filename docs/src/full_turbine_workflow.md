# Full Turbine Workflow

The high-level turbine workflow is the interface used by coupled OWENS simulations. It builds stacked aerodynamic slices once, stores solver state, and advances or re-solves the aero model as the structural state changes.

## Lifecycle

1. `setupTurb(...)` constructs slice geometry, reads airfoil data, initializes environments, and creates unsteady state containers.
2. `steadyTurb(; omega, Vinf)` evaluates the initialized turbine at a steady operating point.
3. `advanceTurb(tnew; ts, azi, verbosity, alwaysrecalc, last_step)` advances the unsteady solution to a new time and returns load histories and state.
4. `AdvanceTurbineInterpolate(...)` interpolates between adjacent azimuthal evaluations for coupled time integration.

These functions are intentionally stateful. A new call to `setupTurb` is required when rotor geometry, airfoils, model options, inflow source, or discretization change.

## Geometry Inputs

`setupTurb` expects blade centerline coordinates by slice, blade count, chord, rotation rate, and inflow. The examples use `shapeX` for radius-like positions and `shapeZ` for vertical station positions. Chord can be scalar or station-dependent depending on the call path.

Airfoils are supplied with `afname`. The Boeing-Vertol reader is used when `DynamicStallModel = "BV"` and the simpler reader is used when dynamic stall is disabled. For single-Reynolds-table reads, the BV reader uses the stall-angle, zero-lift-angle, and thickness metadata in the polar header. Polar files with a `Cm25` column propagate that coefficient into the returned `cm_af` diagnostics and the `M25` distributed pitching-moment array.

`setupTurb` passes `speed_of_sound` into each slice `Environment`, where DMS and AC use it to form the Mach number passed to airfoil callbacks. The default is `343.0` m/s for legacy air cases; water and temperature-specific validation cases should set it explicitly.

For helical blades, pass the optional `bld_y` centerline coordinates alongside `bld_x` and `bld_z`. `setupTurb` converts the local `atan(bld_y, bld_x)` phase into integer azimuth-bin offsets for each slice, with the first node anchored at the blade-root connection. `advanceTurb` applies that same helical offset when it selects blade loads and the scalar per-slice diagnostics (`cl`, `cd_af`, `cm_af`, and `Re`) so those diagnostics follow the shifted blade section rather than the unshifted revolution step.

## Output Shape

The high-level functions return power, torque, per-slice distributed loads, azimuthal angles, local velocities, angles of attack, airfoil coefficients, Reynolds numbers, airfoil pitching moments, and model-specific auxiliary loads. The exact tuple differs between steady and unsteady workflows, so tests should destructure by the current function contract instead of assuming a shared tuple layout.

For unsteady CP, torque, or load-history averages, use `wholeRevolutionIndexRange` or `wholeRevolutionMean` on the continuously increasing azimuth history. These helpers select complete terminal revolutions and exclude the repeated terminal phase so validation metrics are not biased by partial cycles or duplicated boundary samples.

## Coupling Notes

When coupled to OWENS structural dynamics:

- update `setupTurb` only for topology or discretization changes;
- use `deformTurb` or the coupled wrapper to update blade deflection state;
- keep aero and structural sign conventions visible in the test name and fixture;
- verify `M25` separately when pitching-moment polars are used, because it maps to the structural DOF-4 distributed moment channel;
- treat `towerShadowVelocity` as a caller-side inflow primitive until tower-shadow feedback is validated in the solver loop;
- map lumped ancillary loads such as `jointDragForce` explicitly in the caller until a validated solver-integrated drag option exists;
- map lifting-member loads from `liftingStrutForce` explicitly in the caller until strut geometry, projected lift directions, airfoil tables, and structural load transfer are validated end to end;
- pin force components separately instead of testing only aggregate power.

The examples in `test/full_turb.jl`, `test/full_turb_undersampling.jl`, and `examples/RM2/RM2_medium.jl` are the most useful starting points for complete turbine cases.
