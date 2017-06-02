using Documenter, PDMP

try
    # refreshes all the examples
    include("readexamples.jl")
catch e
    # probably on travis CI
    println(e)
end

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
    assets = String[
                "assets/partial.css",
                "assets/BPS.svg",
                "asserts/truncatedgaussian.png"],
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDMP.jl.git",
    deps   = Deps.pip("mkdocs", "pygments"),
    target = "build",
    julia  = "0.5",
    deps = nothing,
    make = nothing
)
