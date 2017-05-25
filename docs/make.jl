using Documenter, PDMP

makedocs(
    modules = [PDMP],
    format = :html,
    sitename = "PDMP.jl",
    authors = "Thibaut Lienart",
    pages = Any[
        "Introduction" => "index.md",
        "About PDMP" => "aboutpdmp.md",
        "Examples" => Any[
            "Global BPS" => "examples/ex_gbps1.md",
            "Local BPS" => "examples/ex_lbps1.md"
        ],
        "Technical Documentation" => Any[
            "Code structure" => "techdoc/structure.md",
            "Types" => "techdoc/types.md"
        ],
        "Contributing" => Any[
            "Examples" => "contributing/addingexample.md"
        ]
    ],
    assets = String["assets/partial.css", "assets/BPS.svg"],
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    deps   = Deps.pip("mkdocs", "pygments"),
    target = "build",
    julia  = "0.5",
    deps = nothing,
    make = nothing
)
