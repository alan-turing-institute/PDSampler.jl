using Documenter, PDMP

makedocs()

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    branch = "gh-pages",
    deps   = Deps.pip("mkdocs","python-markdown-math"),
    julia  = "0.5"
    )
