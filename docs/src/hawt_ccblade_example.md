# Rigid HAWT CCBlade Example

`examples/HAWT/rigid_ccblade_hawt.jl` is the current plot-free horizontal-axis
rotor example. It uses the exported CCBlade adapter as the steady BEM backend
and then advances one Oye dynamic-inflow state update from the steady axial
induction.

The example intentionally stays rigid and aerodynamic only. It does not map
loads into OWENSFEA, rotate a stationary tower in the global frame, or compare
against OpenFAST yet. Those are the next validation steps before this path
should be treated as an aeroelastic HAWT solver.

Run it from the package root with:

```julia
julia --project=. examples/HAWT/rigid_ccblade_hawt.jl
```

The returned metrics are also pinned in `test/hawt_example_tests.jl`:

- `CP = 0.1774837007705498`
- `CT = 0.22007046759243654`
- `CQ = 0.10141925758317133`
- `thrust = 1327.9868847915095 N`
- `torque = 4284.01010759976 N*m`
- `power = 8568.02021519952 W`

The default HAWT solve uses CCBlade's Prandtl tip-plus-hub correction. The API
also exposes `tip_correction = nothing` and `tip_correction =
CCBlade.PrandtlTip()` for comparison studies, and the tests pin the returned
loss factors against OWENSAero's `prandtlTipLossFactor` helper.

The AeroDyn reader helpers can import the blade geometry and AirfoilInfo polar
tables used by OpenFAST verification cases:

```julia
primary = readAeroDynPrimaryFile("ad_primary.dat")
blade = readAeroDynBladeFile("NRELOffshrBsline5MW_AeroDyn_blade.dat")
polar = readAeroDynAirfoilInfo("Airfoils/DU40_A17.dat")
af = aeroDynAirfoilFunction(polar)
```

Those helpers are input-normalization tools, not a completed OpenFAST
validation by themselves. The comparison still needs an explicit root-station
policy, wind-source selection, drag-in-induction setting, and torque/power sign
convention.

The Oye state update uses the steady CCBlade axial induction as the
quasi-steady input. This keeps the future dynamic-inflow coupling explicit:
the HAWT solver must own induction state, time constants, and load mapping
rather than reusing the DMS scalar wake-speed filter.
