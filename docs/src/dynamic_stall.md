# Dynamic Stall

OWENSAero exposes the Boeing-Vertol dynamic-stall model through
`Boeing_Vertol` and through the high-level turbine workflow when
`DynamicStallModel = "BV"`.

The CI regression test in `test/dyn_stall_tests.jl` runs the NACA 0012
pitching case derived from Sadr et al. figure data. It pins the full
angle-of-attack, lift, drag, and pitching-moment histories against the
checked-in HDF5 fixture and also checks representative extrema and endpoint
values. The test intentionally avoids GUI plotting and optional plotting
packages.

For a compact runnable example without plotting dependencies:

```bash
julia --project=. examples/dynamic_stall/BoeingVertol_NACA0012.jl
```

That example uses an analytic static airfoil callback to show the state update
pattern. For validation work, use measured or published static polar data and
pin the resulting histories or scalar metrics in tests.
