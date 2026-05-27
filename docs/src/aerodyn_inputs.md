# AeroDyn Input Readers

OWENSAero includes lightweight readers for AeroDyn v15 blade-definition files
and AirfoilInfo v1 polar files, plus the primary AeroDyn option file. These
helpers support the HAWT validation path: they make the OpenFAST/AeroDyn BEM
settings, blade geometry, and polar inputs visible to the CCBlade adapter
without routing through the native AeroDyn wrapper.

```julia
primary = readAeroDynPrimaryFile("ad_primary.dat")
blade = readAeroDynBladeFile("NRELOffshrBsline5MW_AeroDyn_blade.dat")
polar = readAeroDynAirfoilInfo("Airfoils/DU40_A17.dat")
airfoil = aeroDynAirfoilFunction(polar)

sections = ccbladeHAWTSections(
    blade.span[2:end] .+ hub_radius,
    blade.chord[2:end],
    blade.twist_rad[2:end],
    airfoil,
)
```

For a complete steady CCBlade solve from the same files:

```julia
result = ccbladeHAWTSolveFromAeroDyn(
    "ad_primary.dat",
    rpm * 2pi / 60,
    inflow_speed,
    rho;
    root_station_policy = :drop_zero_span,
    hub_radius = hub_radius,
)
```

For AeroDyn-driver validation fixtures, prefer reading the driver operating
point instead of copying density, wind speed, rotor speed, and pitch into a
test script:

```julia
driver = readAeroDynDriverFile("HAWT_standalone_test.dvr")
result = ccbladeHAWTSolveFromAeroDynDriver(
    "HAWT_standalone_test.dvr";
    root_station_policy = :drop_zero_span,
    hub_radius = hub_radius,
)
```

Use `aeroDynHAWTCCBladeInputs(...)` when you need to inspect the parsed
geometry, selected airfoil files, station indices, or comparison notes before
running the solve.

`readAeroDynPrimaryFile` extracts the steady BEM and airfoil-table controls
needed for comparison studies: `WakeMod`, `AFAeroMod`, tip and hub loss flags,
tangential induction, induction-drag flags, table column mapping, airfoil-file
paths, blade-file paths, and `UseBlCm`. Paths are preserved as written; callers
should decide whether they are relative to the driver, primary file, or current
working directory.

`readAeroDynDriverFile` extracts the steady HAWT operating point from
`AnalysisType = 1` AeroDyn-driver files. It supports direct steady wind
(`CompInflow = 0`) and steady InflowWind files (`CompInflow = 1`,
`WindType = 1`), and returns the resolved AeroDyn/InflowWind filenames,
density, dynamic viscosity, sound speed, rotor speed, uniform blade pitch, and
driver hub-radius metadata. `BasicHAWTFormat=true` drivers are supported; their
`HubRad` value is used as the physical blade-root radius. It intentionally
rejects non-steady inflow and nonuniform blade pitch because those need a
time-marching validation harness rather than a single rigid CCBlade solve.

`readAeroDynBladeFile` consumes exactly `NumBlNds` blade rows. This matters for
the checked NREL 5 MW fixture, which carries a historical note and an extra
numeric row after the declared blade-node block. The reader validates finite
geometry, strictly increasing span, positive chord, and positive integer
airfoil IDs, then returns both degree and radian forms for twist and curve
angle.

`readAeroDynAirfoilInfo` retains the static `Cl`, `Cd`, and `Cm` columns plus
the unsteady-aero metadata needed for later dynamic-stall validation. The
callable produced by `aeroDynAirfoilFunction` uses radians for angle of attack
and returns `(cl, cd)` by default or `(cl, cd, cm)` with `return_cm=true`.
Queries outside the tabulated angle range fail instead of silently
extrapolating.

For HAWT comparisons, keep the frame conventions explicit:

- AeroDyn blade twist is stored in degrees; the CCBlade adapter expects radians.
- AeroDyn blade files store `BlSpn` from the blade root, while CCBlade sections
  use rotor radius from the center of rotation. `aeroDynHAWTCCBladeInputs`
  therefore adds `hub_radius` to each retained span station.
- AeroDyn blade files can include a root station at zero span. The validation
  setup must state whether it drops the root station or applies a documented
  hub/root policy.
- AeroDyn BEM option flags are not one-to-one with CCBlade defaults. In
  particular, drag-in-induction and hub-loss settings need explicit matching
  before comparing coefficients.
- OpenFAST driver files may delegate wind to InflowWind; do not assume the
  driver's `HWindSpeed` is the operating wind speed when `CompInflow = 1`.
- For station outputs in the Basic HAWT validation fixture, CCBlade `cn`/`ct`
  map to AeroDyn `Cx`/`Cy`, and CCBlade `Np`/`Tp` map to AeroDyn `Fxp`/`-Fyp`.
  Root and exact-tip stations are treated separately because hub/tip-loss
  behavior is singular there.
