using Documenter, Literate, VAWTAero

# Build documentation
makedocs(;
    modules = [VAWTAero],
    pages = [
        "Home" => "index.md",
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "VAWTAero.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
)