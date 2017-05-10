using Documenter, PDMP

makedocs()

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    target = "build",
    branch = "gh-pages",
    julia  = "0.5",
    osname = "linux",
    deps   = Deps.pip("pygments", "mkdocs", "python-markdown-math")
)
