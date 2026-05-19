# Dynamic Stall

OWENSAero exposes the Boeing-Vertol dynamic-stall model through
`Boeing_Vertol` and through the high-level turbine workflow when
`DynamicStallModel = "BV"`.

When `readaerodyn_BV` or `setupTurb(...; DynamicStallModel = "BV")` reads an
OWENSAero polar file, the Boeing-Vertol positive stall angle, negative stall
angle, zero-lift angle, and thickness-to-chord ratio come from the first
Reynolds-number table in that file. This keeps dynamic-stall state transitions
tied to the checked-in polar metadata rather than hidden reader defaults.

The CI regression test in `test/dyn_stall_tests.jl` runs the NACA 0012
pitching case derived from Sadr et al. figure data. It pins the full
angle-of-attack, lift, drag, and pitching-moment histories against the
checked-in HDF5 fixture and also checks representative extrema and endpoint
values. It also computes branch-wise upstroke/downstroke lift and
pitching-moment errors against the digitized Sadr Boeing-Vertol figure data.
Those metrics make dynamic lift and moment drift visible in CI without relying
on GUI plots or optional plotting packages.

The checked-in validation data do not include a digitized dynamic drag figure.
The drag history is therefore regression-pinned against the HDF5 fixture and
representative scalar values, but coefficient tuning or acceptance claims for
dynamic `CD` should wait until measured, published, or legacy-reference dynamic
drag data are added to the test set.

For a compact runnable example without plotting dependencies:

```bash
julia --project=. examples/dynamic_stall/BoeingVertol_NACA0012.jl
```

That example uses an analytic static airfoil callback to show the state update
pattern. For validation work, use measured or published static polar data and
pin the resulting histories or scalar metrics in tests.
