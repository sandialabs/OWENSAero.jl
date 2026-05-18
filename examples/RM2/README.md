# RM2 Example

`RM2_medium.jl` is a standalone OWENSAero stacked-slice RM2 case with the
checked-in RM2 airfoils and experimental CP data. It runs the DMS branch by
default so it stays practical as an example:

```bash
julia --project=. examples/RM2/RM2_medium.jl
```

To include both DMS and AC in the same sweep:

```bash
OWENSAERO_RM2_MODELS=all julia --project=. examples/RM2/RM2_medium.jl
```

`PlotQbladeResults.jl` overlays the checked-in RM2 experimental CP curve with
representative QBlade points and writes the figure under `examples/RM2/figures`,
which is ignored as generated output.
