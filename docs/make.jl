using Documenter, PDMP

makedocs(
    modules = [PDMP],
    format  = :html,
    sitename = "PDMP.jl",
    authors = "Thibaut Lienart",
    pages   = Any[
        "Home"     => "index.md",
        "Examples" => Any[
            "Global BPS 1"  => "examples/bps_mvg_constr.md",
            "Local BPS 1"   => "examples/lbps_gchain.md"
        ],
        "Technical doc" => [
            "Types" => "techdoc/types.md"
        ],
    ],
    doctest = false
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    target = "build",
    branch = "gh-pages",
    julia  = "0.5",
    osname = "linux",
    deps   = Deps.pip("pygments", "mkdocs", "python-markdown-math")
)
