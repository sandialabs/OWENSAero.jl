# Aero Model Audit

This page records the current implementation and verification state for the
major OWENSAero submodels. It is intended to keep model selection, validation
evidence, and open research work aligned.

## Verification Standard

New submodels should include:

- a pure API test with concrete numerical pins and invalid-input coverage;
- a solver-level regression showing how the option changes loads or confirms
  that it is intentionally off by default;
- documentation that states whether the model is coupled into DMS, AC, or only
  available as a caller-side primitive;
- validation metrics against experiment, a trusted code, or an analytic limit
  before changing existing DMS/AC baselines.

## Current Model State

| Model or addition | Current status | Verification in tree | Remaining work |
| --- | --- | --- | --- |
| Double-multiple streamtube (DMS) | Coupled solver path for stacked two-dimensional VAWT slices. | Regression tests pin signed CP, negative RPM behavior, wind-angle invariance, helical index mapping, and added-mass load channels. | Add finite-height correction only behind an explicit option with validation; continue high-TSR/high-solidity diagnostics. |
| Actuator cylinder (AC) | Coupled public path for a single turbine. Internal multi-turbine assembly helpers exist, but the public `AC` path rejects multi-turbine calls. | Regression tests pin single-turbine AC force channels, moment propagation, added mass, and negative RPM behavior. | Restore and validate public multi-turbine AC before documenting it as usable; add gradient and radial-influence validation coverage. |
| Static polar lookup and moment coefficient | Coupled through DMS/AC and propagated to structural moment channels through the OWENS stack. | Tests pin `Cm25` propagation, polar parsing, and OWENS structural mapping. | Unify OWENSAero, AeroDyn, XFOIL, and WindIO polar formats; document and validate extrapolation policy outside tabulated alpha/Re ranges. |
| Boeing-Vertol dynamic stall | Implemented and coupled through the airfoil function path. | Headless tests pin CL/CM branch behavior and validation metrics against digitized Sadr dynamic-stall data. | Add a real dynamic-CD validation dataset before changing drag behavior. |
| Leishman-Beddoes dynamic stall | Not implemented as a production option. `DynamicStallModel = "LB"` fails fast. | Constructor/setup tests pin the fail-fast behavior. | Implement state equations, lifecycle, tests, and validation metrics before exposing it. |
| Finite-blade tip and root loss | `prandtlTipLossFactor` is an exported, tested axial-flow primitive. DMS/AC still use the no-tip-loss stacked-slice baseline. | API tests pin midspan, tip, root, near-zero inflow angle, invalid inputs, and ForwardDiff derivative behavior. | Decide and validate a VAWT finite-height correction before coupling any loss factor into DMS or AC. |
| Added mass, buoyancy, and rotational acceleration | Coupled through DMS/AC options and mapped into OWENS structural load channels. | Zero-aero tests pin geometry factors, force signs, mass channels, and structural mapping. | Add physical validation metrics for force-coupled added mass and buoyancy before changing equations or signs. |
| Joint drag | Exported helper only. It does not modify DMS/AC induction or structural loads unless a caller maps it. | API tests pin vector direction, magnitude, zero cases, and invalid inputs. | Add geometry mapping and structural load-channel integration with validation data. |
| Tower shadow | Exported helper only. It does not modify DMS/AC induction unless a caller applies the returned velocity. | API tests pin off-by-default behavior, downstream gating, lateral deficit, wake expansion, and invalid inputs. | Add a solver option and validate against tower-shadow or wake-deficit data before changing blade loads. |
| Lifting struts | Exported helper only. It does not modify DMS/AC induction or OWENS structural channels unless a caller maps it. | API tests pin 2D/3D lift-drag decomposition, oblique directions, zero cases, and invalid inputs. | Couple strut geometry, aero data, postprocessing, and VTK/stress output after validation. |
| Unsteady wake and RPI | Implemented for single-turbine positive-rotation workflows. Negative mean rotor speed with `RPI=true` fails fast until the clockwise unsteady path is validated. | Examples and tests pin representative unsteady histories, whole-revolution averaging windows, wake-speed floor behavior, positive `tau` validation, and the negative-RPM RPI guard. | Add timestep-sensitivity and aliasing diagnostics; validate clockwise and multi-rotor behavior before expanding support. |
| HAWT Oye dynamic inflow | Exported helper only. It advances reduced and dynamic induction states for future CCBlade-based HAWT work. | API tests pin Oye time constants, state derivatives, exact constant-input steps, equal-time-constant handling, invalid inputs, and a ForwardDiff derivative. | Add the CCBlade adapter, rigid-rotor HAWT example, and OpenFAST/AeroDyn DBEMT validation before coupling HAWT loads into OWENS. |
| InflowWind coupling | Available through high-level setup paths when external TurbSim/OpenFAST inputs are supplied. | Existing fixtures cover required file layout and selected setup behavior. | Harden WindIO airfoil usage and secondary-file provenance in validation examples. |

## Additional Model Work To Triage

The next aero-model issues should be tracked explicitly instead of being folded
into existing solvers without evidence:

- validated VAWT finite-height/tip-loss correction for DMS and AC;
- dynamic-CD validation for Boeing-Vertol and a full Leishman-Beddoes model;
- restored public multi-turbine AC with regression metrics;
- polar-format unification and out-of-range interpolation/extrapolation policy;
- solver-integrated tower shadow, joint drag, and lifting-strut loads;
- timestep/azimuth-convergence diagnostics for unsteady DMS and AC;
- HAWT/axial-flow aero path using CCBlade, Oye dynamic inflow, global
  rotation-frame bookkeeping, and OpenFAST comparison metrics.

## API Reference

The exported primitives on this page, including `prandtlTipLossFactor`, are
also pulled into the generated API reference through the package autodocs.
