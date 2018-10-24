using Documenter, PDSampler

try
    # refreshes all the examples
    include("readexamples.jl")
catch e
    # probably on travis CI
    println(e)
end

makedocs(
    modules = [PDSampler],
    format = :html,
    sitename = "PDSampler.jl",
    authors = "Thibaut Lienart",
    pages = Any[
        "Introduction" => "index.md",
        "About PDSampler"   => "aboutpdsampler.md",
        "Examples" => Any[
            "Global BPS" => "examples/ex_gbps1.md",
            "Local BPS"  => "examples/ex_lbps1.md"
        ],
        "Technical Documentation" => Any[
            "Code structure" => "techdoc/structure.md",
            "Core tools"     => "techdoc/coretools.md",
            "Models"         => "techdoc/models.md",
            "Global sampler" => "techdoc/global.md",
            "Local sampler"  => "techdoc/local.md",
#            "Types"          => "techdoc/types.md",
        ],
        "Contributing" => Any[
            "New example" => "contributing/addingexample.md"
            "New feature" => "contributing/addingfeature.md"
        ]
    ],
    assets = String[
                "assets/partial.css",
                "assets/BPS.svg",
                "asserts/truncatedgaussian.png"],
)

deploydocs(
    repo   = "github.com/alan-turing-institute/PDSampler.jl.git",
    deps   = Deps.pip("mkdocs", "pygments"),
    target = "build",
    julia  = "0.5",
    make = nothing
)
