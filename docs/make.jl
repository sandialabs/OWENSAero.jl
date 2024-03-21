using Documenter, Literate, OWENSAero

# Build documentation
makedocs(;
    modules = [OWENSAero],
    pages = [
        "Home" => "index.md",
        "API Reference" => joinpath("reference", "reference.md"),
    ],
    sitename = "OWENSAero.jl",
    authors = "Kevin R. Moore <kevmoor@sandia.gov>",
)