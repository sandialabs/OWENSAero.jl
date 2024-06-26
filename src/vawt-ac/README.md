# VAWT Multiple Actuator Cylinder Model

[![DOI](https://zenodo.org/badge/51314790.svg)](https://zenodo.org/badge/latestdoi/51314790)

Author: Andrew Ning

Major invasive changes for OWENS code.  Author: Kevin Moore
- Add Dynamic Stall Model from CACTUS code, interface to airfoil reading methods
- RPI method (see paper)
- Corrections for curved and deforming blades, turbulence, and non-constant rotational speed.
- Modifications to propagate automatic gradients
- Modification to handle non-x direction wind

A modified version of actuator cylinder theory to account for multiple vertical axis wind turbines.  Can run with just a single turbine.  See examples in example.jl.

Written in Julia, but an older Python version of the single actuator cylinder theory is also available in a branch.

See paper for theory details: Ning, A., "Actuator Cylinder Theory for Multiple Vertical Axis Wind Turbines," Wind Energy Science, Nov 2016, (accepted). doi:10.5194/wes-2016-19
