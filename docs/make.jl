using Documenter, Literate, OWENSAero

# Build documentation
makedocs(;
    modules = [OWENSAero],
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Model Selection" => "model_selection.md",
        "Full Turbine Workflow" => "full_turbine_workflow.md",
        "Dynamic Stall" => "dynamic_stall.md",
        "RM2 Example" => "rm2_example.md",
        "Frames, Units, and Outputs" => joinpath("theory", "frames_units.md"),
        "Validation and Testing" => "validation.md",
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSAero.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
    remotes = nothing,
)

deploydocs(
    repo = "github.com/sandialabs/OWENSAero.jl.git",
)
