# RM2 Example

The RM2 example under `examples/RM2` is the standalone aerodynamic companion to
the broader OWENS RM2 workflow. It uses the checked-in RM2 airfoils and
experimental CP data and writes generated figures under an ignored `figures`
directory.

Run the default DMS sweep with:

```bash
julia --project=. examples/RM2/RM2_medium.jl
```

Include both DMS and AC by setting:

```bash
OWENSAERO_RM2_MODELS=all julia --project=. examples/RM2/RM2_medium.jl
```

The separate `PlotQbladeResults.jl` script overlays representative QBlade CP
points with the checked-in experimental RM2 curve. These scripts are examples,
not CI validation. A CI-quality RM2 validation should pin scalar error metrics
against the experimental data for a specific model configuration and runtime
budget.
