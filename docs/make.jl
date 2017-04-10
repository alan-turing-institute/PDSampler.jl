using Documenter, PDMP

makedocs(
    modules  = [PDMP],
    sitename = "PDMP.jl",
    authors  = "Thibaut Lienart",
    pages    = Any[
        "Home" => "index.md"
    ]
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    branch = "gh-pages",
    deps   = Deps.pip("mkdocs","python-markdown-math"),
    julia  = "0.5"
    )
