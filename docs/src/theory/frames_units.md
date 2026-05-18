# Frames, Units, and Outputs

OWENSAero uses SI units. Angles are radians unless a field or helper name explicitly contains `_d`.

## Coordinate Conventions

| Quantity | Convention |
| --- | --- |
| `R`, `r` | Turbine radius and local radius in meters. |
| `z` or `shapeZ` | Vertical slice or station location in meters. |
| `theta` / azimuth arrays | Azimuthal blade position in radians. |
| `omega` | Rotor angular velocity in radians per second. |
| `windangle` | Mean inflow angle in radians. |
| `V_x`, `V_y`, `V_z` | Velocity components in meters per second. |

For simple runs, `V_x` is the incoming freestream velocity. `V_y` and `V_z` are zero unless cross-flow, vertical flow, or coupled motion terms are active.

## Load Components

The direct slice path returns `Rp`, `Tp`, and `Zp`:

| Output | Meaning | Typical unit |
| --- | --- | --- |
| `Rp` | radial or normal force per height | N/m |
| `Tp` | tangential force per height | N/m |
| `Zp` | vertical force per height | N/m |
| `M25` | quarter-chord airfoil pitching moment per span, computed as `Cm25 * 0.5 * rho * Vloc^2 * chord^2` | N m/m |
| `Q` | torque | N m |
| `CP`, `CD`, `CT`, `cl`, `cd_af`, `cm_af` | nondimensional performance and force coefficients | dimensionless |

The coupled turbine path maps those sectional quantities to blade and global loads before passing them to OWENS. Any new coupled validation should assert both local aerodynamic loads and transformed structural loads to catch sign or frame mistakes.

## Airfoil Data

Airfoil reader functions expect angle-of-attack and coefficient tables compatible with the selected model. Tables may include a fourth `Cm25` column. The readers keep the older `cl, cd = af(alpha, Re, Mach)` call pattern intact, and `af(alpha, Re, Mach; return_cm=true)` returns `cl`, `cd`, and `Cm25`.

DMS and AC convert `Cm25` to a distributed section moment with `M25 = Cm25 * q * chord^2`, where `q = 0.5 * rho * Vloc^2`. That moment is returned with the aerodynamic loads for coupled OWENS runs. If a user-supplied airfoil function only returns `cl, cd`, OWENSAero uses `Cm25 = 0.0`.

The dynamic-stall path carries additional state for prior angle of attack and Boeing-Vertol flags. Because these state variables affect transient loads, unsteady tests should pin the full output vector for a named input file and time history.

## Automatic Differentiation

The tests include gradient checks for representative steady DMS paths. AC gradient coverage is still incomplete, so gradient-sensitive studies should avoid undocumented AC branches and dynamic-stall branches unless the selected branch has explicit derivative tests. Keep airfoil interpolation inputs within tabulated ranges; extrapolation paths are intentionally guarded because undefined spline behavior can corrupt gradients.
