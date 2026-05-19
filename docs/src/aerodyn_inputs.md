# AeroDyn Input Readers

OWENSAero includes lightweight readers for AeroDyn v15 blade-definition files
and AirfoilInfo v1 polar files. These helpers support the HAWT validation path:
they make the OpenFAST/AeroDyn blade and polar inputs visible to the CCBlade
adapter without routing through the native AeroDyn wrapper.

```julia
blade = readAeroDynBladeFile("NRELOffshrBsline5MW_AeroDyn_blade.dat")
polar = readAeroDynAirfoilInfo("Airfoils/DU40_A17.dat")
airfoil = aeroDynAirfoilFunction(polar)

sections = ccbladeHAWTSections(
    blade.span[2:end],
    blade.chord[2:end],
    blade.twist_rad[2:end],
    airfoil,
)
```

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
- AeroDyn blade files can include a root station at zero span, while CCBlade
  sections require positive radial positions. The validation setup must state
  whether it drops the root station or applies a documented hub/root policy.
- OpenFAST driver files may delegate wind to InflowWind; do not assume the
  driver's `HWindSpeed` is the operating wind speed when `CompInflow = 1`.
- Existing wrapper fixtures can report torque and power with the OpenFAST
  output-channel sign convention. Native CCBlade comparisons should normalize
  the torque/power frame before computing validation metrics.
