# Dynamic Stall

OWENSAero exposes the Boeing-Vertol and Leishman-Beddoes dynamic-stall models
through `Boeing_Vertol`, `Leishman_Beddoes`, and through the high-level turbine
workflow when `DynamicStallModel = "BV"` or `DynamicStallModel = "LB"`.

When `readaerodyn_BV` or `setupTurb(...; DynamicStallModel = "BV")` reads an
OWENSAero polar file, the Boeing-Vertol positive stall angle, negative stall
angle, zero-lift angle, and thickness-to-chord ratio come from the first
Reynolds-number table in that file. This keeps dynamic-stall state transitions
tied to the checked-in polar metadata rather than hidden reader defaults.

When `DynamicStallModel = "LB"`, `readaerodyn_BV_NEW` uses the same static
polar table and reads the Leishman-Beddoes lift-curve slope, positive critical
lift coefficient, and negative critical lift coefficient from the airfoil file.
The model stores lagged pressure, trailing-edge separation, leading-edge vortex,
and leading-edge separation-state histories in `LeishmanBeddoesState`. During
DMS and actuator-cylinder residual solves, those histories are read but not
mutated; they are updated only for accepted section evaluations. This preserves
the same solver-state bookkeeping rule used by the Boeing-Vertol path.

The CI regression test in `test/dyn_stall_tests.jl` runs the NACA 0012
pitching case derived from Sadr et al. figure data. It pins the full
angle-of-attack, lift, drag, and pitching-moment histories against the
checked-in HDF5 fixture and also checks representative extrema and endpoint
values. It also computes branch-wise upstroke/downstroke lift and
pitching-moment errors against the digitized Sadr Boeing-Vertol figure data.
Those metrics make dynamic lift and moment drift visible in CI without relying
on GUI plots or optional plotting packages.

The Leishman-Beddoes regression test uses the same NACA 0012 pitching case and
the DynStall.jl-style critical-lift inputs. It pins the full angle-of-attack,
lift, drag, pitching-moment, leading-edge separation-state, and leading-edge
vortex-age histories against a checked-in HDF5 fixture, then reports branch-wise
errors against the digitized Sadr Leishman-Beddoes lift and moment curves.

The checked-in validation data do not include a digitized dynamic drag figure.
The drag history is therefore regression-pinned against the HDF5 fixture and
representative scalar values, but coefficient tuning or acceptance claims for
dynamic `CD` should wait until measured, published, or legacy-reference dynamic
drag data are added to the test set.

For a compact runnable example without plotting dependencies:

```bash
julia --project=. examples/dynamic_stall/BoeingVertol_NACA0012.jl
julia --project=. examples/dynamic_stall/LeishmanBeddoes_NACA0012.jl
```

That example uses an analytic static airfoil callback to show the state update
pattern. For validation work, use measured or published static polar data and
pin the resulting histories or scalar metrics in tests.
