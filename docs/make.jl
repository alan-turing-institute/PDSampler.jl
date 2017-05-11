using Documenter, PDMP

listofexamples = filter(e->e[end-2:end]==".md",readdir("src/examples/"))

makedocs(
    modules = [PDMP],
    format = :html,
    sitename = "PDMP.jl",
    authors = "Thibaut Lienart",
    pages = Any[
        "Introduction" => "index.md",
        "Examples" => Any[
            "examples/ex_gbps1.md",
            "examples/ex_lbps1.md"
        ],
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
