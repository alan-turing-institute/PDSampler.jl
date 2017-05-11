using Documenter, PDMP

include("readexamples.jl")

makedocs(
    modules = [PDMP],
    format = :html,
    sitename = "PDMP.jl",
    authors = "Thibaut Lienart",
    pages = Any[
        "Introduction" => "index.md",
        "Examples" => listofexamples,
        "Technical Documentation" => Any[
            "Types" => "techdoc/types.md"
        ]
    ]
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    target = "build",
    julia  = "0.5",
    deps = nothing,
    make = nothing
)
