using Documenter, PDMP

makedocs(
    modules  = [PDMP],
    doctest  = false,
    sitename = "PDMP.jl",
    authors  = "Thibaut Lienart"
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    branch = "gh-pages",
    deps   = Deps.pip("mkdocs","python-markdown-math"),
    julia  = "0.5"
    )
