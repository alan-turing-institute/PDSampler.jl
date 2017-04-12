__precompile__(false) # TODO change to true when more stable

# =================================================================
# TODO
# at the moment pretty much all functions are listed.
# this is useful to keep track of what's documented/tested
# however in practice, it's not necessary to export all methods
# =================================================================

module PDMP

if Sys.WORD_SIZE == 32
    typealias Int Int32
    typealias Float Float32
elseif Sys.WORD_SIZE == 64
    typealias Int Int64
    typealias Float Float64
else
    error("Unknown architecture...")
end

export
    Int,
    Float,
    # MODELS
    loglik,
    gradloglik,
    gradloglik_cv

### source files (keep the order)

include("models/mvgaussian.jl")
include("models/logreg.jl")
include("models/pmf.jl")

include("geometry.jl")
include("ippsampler.jl")
include("path.jl")
include("kernels.jl")
include("simulate.jl")

# include("local/factor.jl")
include("local/event.jl")
include("local/factorgraph.jl")
include("local/simulate.jl")

end # module
