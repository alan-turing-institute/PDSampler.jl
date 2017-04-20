__precompile__(true) # TODO set to true in releases

# =================================================================
# TODO
# at the moment pretty much all functions are listed.
# this is useful to keep track of what's documented/tested
# however in practice, it's not necessary to export all methods
# =================================================================

module PDMP

const Int   = Int64
const Float = Float64

export
    Int,
    Float,
    # MODELS
    loglik,
    gradloglik,
    gradloglik_cv

using Polynomials:
        Poly,
        roots,
        polyint,
        polyval

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
