using Documenter, Literate, OWENSAero

# Build documentation
makedocs(;
    modules = [OWENSAero],
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Model Selection" => "model_selection.md",
        "Aero Model Audit" => "model_audit.md",
        "Full Turbine Workflow" => "full_turbine_workflow.md",
        "Dynamic Stall" => "dynamic_stall.md",
        "RM2 Example" => "rm2_example.md",
        "AeroDyn Input Readers" => "aerodyn_inputs.md",
        "Rigid HAWT CCBlade Example" => "hawt_ccblade_example.md",
        "Frames, Units, and Outputs" => joinpath("theory", "frames_units.md"),
        "Validation and Testing" => "validation.md",
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSAero.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
    remotes = nothing,
    format = Documenter.HTML(
        repolink = "https://github.com/sandialabs/OWENSAero.jl",
        edit_link = "master",
    ),
)

if get(ENV, "CI", "false") == "true"
    deploydocs(repo = "github.com/sandialabs/OWENSAero.jl.git", devbranch = "master")
end
