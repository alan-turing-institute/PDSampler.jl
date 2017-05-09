using Documenter, PDMP

makedocs(
    modules  = [PDMP],
    doctest  = false,
    clean    = true,
    format   = :html,
    sitename = "PDMP.jl",
    authors  = "Thibaut Lienart",
    pages    = Any[
        "Home"      => "index.md",
        "Examples"  => Any[
            "examples/bps_mvg_constr.md",
            "examples/lbps_gchain.md"
        ],
        "Technical Documentation" => Any[
            "techdoc/types.md"
        ]
    ]
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    target = "build",
    osname = "linux",
    branch = "gh-pages",
    julia  = "0.5"
    deps   = Deps.pip("mkdocs","python-markdown-math")
    )
