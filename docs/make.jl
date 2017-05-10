using Documenter, PDMP

makedocs(
    modules = [PDMP],
    format  = :html,
    sitename = "PDMP.jl",
    pages   = [
        "Home" => "index.md",
        "Examples" => [
            "Global BPS 1" => "examples/bps_mvg_constr.md",
            "Local BPS 1" => "examples/lbps_gchain.md"
        ],
        "Technical doc" => [
            "Types" => "techdoc/types.md"
        ]
    ]
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    target = "build",
    make   = "nothing",
    julia  = "0.5",
    osname = "linux",
    deps   = Deps.pip("mkdocs", "python-markdown-math"))

# makedocs(
#     modules  = [PDMP],
#     doctest  = false,
#     clean    = true,
#     format   = :html,
#     sitename = "PDMP.jl",
#     authors  = "Thibaut Lienart"
# )
#
# deploydocs(
#     repo   = "github.com/alan-turing-institute/PDMP.jl.git",
#     target = "build",
#     osname = "linux",
#     branch = "gh-pages",
#     julia  = "0.5",
#     deps   = Deps.pip("mkdocs","python-markdown-math")
#     )
