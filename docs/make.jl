using Documenter, PDMP

makedocs(
    modules  = [PDMP],
    doctest  = false,
    clean    = true,
    format   = :html,
    sitename = "PDMP.jl",
    authors  = "Thibaut Lienart"
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    target = "build",
    osname = "linux",
    branch = "gh-pages",
    julia  = "0.5",
    deps   = Deps.pip("mkdocs","python-markdown-math")
    )
