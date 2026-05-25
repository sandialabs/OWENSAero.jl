# Quickstart

This page shows the shortest path from package installation to a pinned aerodynamic result. The examples use SI units throughout.

## Install

For released use, install the registered package or the Sandia repository:

```julia
using Pkg
Pkg.add(PackageSpec(url = "https://github.com/sandialabs/OWENSAero.jl.git"))
```

For toolkit development from sibling checkouts:

```julia
using Pkg
Pkg.develop(path = "../OWENSOpenFASTWrappers.jl")
Pkg.develop(path = "../OWENSAero.jl")
```

The docs project uses local package sources and should be built with Julia 1.11
or newer.

## Single-Slice Solve

The direct slice API is the best starting point for understanding the aerodynamic state. It mirrors the coverage tests in `test/simple_example.jl`.

```julia
using OWENSAero

ntheta = 30
theta = range(0.0, 2pi, length = ntheta + 1)[1:end-1]
R = 1.0
r = fill(R, ntheta)
chord = fill(0.08, ntheta)
twist = zeros(ntheta)
delta = zeros(ntheta)
omega = fill(12.0, ntheta)
B = 3
rho = 1.225
mu = 1.81e-5
Vinf = fill(8.0, ntheta)

af(alpha, Re, M) = (2pi * alpha, 0.01, 0.0)
turbine = Turbine(R, r, chord, twist, delta, omega, B, af, ntheta, false)
env = Environment(rho, mu, Vinf, "none", "DMS", zeros(2ntheta))

result = OWENSAero.steady(turbine, env)
```

For both DMS and AC, `result` contains power/thrust/torque quantities, distributed normal/tangential/vertical loads, local velocity, induction, angle-of-attack, coefficient, Reynolds-number, and added-mass/buoyancy output slots. The exact tuple length is model-specific and is pinned in `test/simple_example.jl`; destructure against the tested contract for the selected model.

## Full-Turbine Solve

Use the stateful high-level workflow when the rotor is represented by multiple vertical slices or when the aero model is coupled to OWENS structural motion.

```julia
using OWENSAero

B = 3
shapeX = fill(1.0, 30)
shapeZ = collect(range(0.0, 10.0, length = 30))
chord = 0.08
omega = 12.0
Vinf = 8.0

OWENSAero.setupTurb(shapeX, shapeZ, B, chord, omega, Vinf;
    Nslices = 1,
    ntheta = 30,
    DynamicStallModel = "none",
    AeroModel = "DMS",
    speed_of_sound = 343.0,
)

result = OWENSAero.steadyTurb(; omega, Vinf)
```

The high-level functions return a large tuple of distributed loads, coefficients, velocities, moments, power, torque, and auxiliary load terms. They also store global model state for speed and OpenFAST-style coupling. Re-run `setupTurb` whenever geometry, model choice, airfoil source, or slice count changes.

## First Checks

After installation, run:

```julia
using Pkg
Pkg.test("OWENSAero")
```

The current regression suite includes DMS and AC slice checks, wind-direction cases, dynamic-stall cases, InflowWind-driven unsteady cases, DMS gradient checks, and full-turbine smoke tests.
