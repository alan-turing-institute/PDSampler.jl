using Documenter, PDMP

makedocs(
    modules = [PDMP],
    clean = false,
    format = :html,
    sitename = "PDMP.jl",
    authors = "Thibaut Lienart",
    pages = Any[
        "Introduction" => "index.md",
        "Examples" => Any[
            "Global BPS 1" => "examples/bps_mvg_constr.md",
            "Local BPS 1" => "examples/lbps_gchain.md"
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
#    branch = "gh-pages",
#    osname = "linux",
#    deps   = Deps.pip("pygments", "mkdocs", "python-markdown-math")
