__precompile__(true)

module PDSampler

using Compat
using Klara.ess
using ApproxFun
using Polynomials:
        Poly,
        roots,
        polyint,
        polyval
using DataStructures:
        PriorityQueue,
        enqueue!,
        dequeue!,
        peek
using ProgressMeter
using Distributions: Beta


const Int   = Int64
const Float = Float64

export
    Int,
    Float,
    # MODELS
    loglik,
    gradloglik,
    gradloglik_cv

### source files (keep the order)

include("models/mvgaussian.jl")
include("models/tmvgaussian.jl")
include("models/logreg.jl")
include("models/pmf.jl")

include("geometry.jl")
include("ippsampler.jl")
include("path.jl")
include("kernels.jl")
include("simulate.jl")

include("local/event.jl")
include("local/factorgraph.jl")
include("local/simulate.jl")

end # module
