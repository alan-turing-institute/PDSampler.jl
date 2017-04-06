__precompile__(false) # TODO change to true when more stable

# =================================================================
# TODO
# at the moment pretty much all functions are listed.
# this is useful to keep track of what's documented/tested
# however in practice, it's not necessary to export all methods
# =================================================================

module PDMP

using Distributions:
        Beta
using Base.Collections:
        PriorityQueue,
        enqueue!,
        dequeue!,
        peek
using Polynomials:
        Poly,
        roots,
        polyint,
        polyval

typealias Int Int64
typealias Float Float64

export
    # IPPSAMPLER.jl
    NextEvent,          # not tested, not doc'd FIXME
    LinearBound,        # not tested, not doc'd FIXME
    # -- functions
    nextevent_bps,      # not tested, not doc'd FIXME
    nextevent_zz,       # not tested, not doc'd FIXME

    # GEOMETRY.jl
    Unconstrained,      # trivial, not doc'd FIXME
    Polygonal,          # trivial, not doc'd FIXME
    # -- functions
    nextboundary,       # not tested, not doc'd FIXME

    # PATH.jl
    AllowedTimeType,    # trivial [has doc]
    Path,               # trivial, not doc'd FIXME
    # -- functions
    samplepath,         # [has test], [has doc]

    # KERNELS.jl
    # -- functions
    reflect_bps!,       # [has test], [has doc]
    reflect_zz!,        # [has doc]
    refresh_global!,    # not tested, not doc'd FIXME
    refresh_restricted!,# not tested, not doc'd FIXME
    refresh_partial!,   # not tested, not doc'd FIXME

    # SIMULATE.jl
    Simulation,         # not tested, not doc'd FIXME
    simulate,           # not tested, not doc'd FIXME

    # MODELS
    loglik,             # [has test], [has doc]
    gradloglik,         # [has test], [has doc]
    gradloglik_cv,

    # MODELS/MVGAUSSIAN.jl
    MvGaussian,         # abstract type
    MvGaussianStandard, # [has test], [has doc]
    MvDiagonalGaussian, # [has test], [has doc]
    MvGaussianCanon,    # [has test], [has doc]
    MvGaussianNatural,  # [has test], [has doc]
    # -- functions
    mvg_mu,             # [has test], [has doc]
    mvg_precmu,         # [has test], [has doc]
    mvg_precmult,       # [has test], [has doc]

    # MODELS/LOGREG.jl
    LogReg,             # [has test], [has doc]
    # -- functions
    logistic,           # [has test], [has doc]
    loglogistic,        # [has test], [has doc]
    gradloglogistic,    # [has test], [has doc]

    # MODELS/PMF.jl
    PMFGaussian,        # not tested, not doc'd FIXME
    # -- functions
    pmf_base,           # not tested, not doc'd FIXME
    pmf_caseA,          # not tested, not doc'd FIXME
    pmf_caseB,          # not tested, not doc'd FIXME
    pmf_caseC,          # not tested, not doc'd FIXME
    pmf_caseD,          # not tested, not doc'd FIXME

    # LOC_SIMULATE
    LocalSimulation,    # not tested, not doc'd FIXME
    # -- function
    simulate,           # not tested, not doc'd FIXME
    ls_reshape,         # [has test], [has doc]
    ls_retrieve,        # [has 1D test], [has doc]
    ls_saveupdate!,     # [has 1D test], [has doc]
    ls_updatepq!,       # [has test], [has doc]
    ls_random,          # [has test], [has doc]
    ls_refreshment,     # not tested, [has partial doc]

    # LOC_FACTORGRAPH
    Factor,             # [has test], [has doc]
    FactorObject,       # not tested, not doc'd FIXME
    FactorGraphStruct,  # [has test], [has doc]
    FactorGraph,        # [has test], [has doc]
    # -- functions
    assocvariables,     # [has test], [has part doc]
    assocfactors,       # [has test], [has part doc]
    linkedfactors,      # [has test], [has doc]
    chainstruct,        # [has test], not doc'd FIXME
    chaingraph,         # [has test], not doc'd FIXME

    # LOC_EVENT
    AllowedVarType,     # [has test], [has doc]
    Event,              # [has test], [has doc]
    EventList,          # [has test], [has doc]
    AllEventList,       # [has test], [has doc]
    # -- functions
    getevent,           # [has test], [has doc]
    getlastevent,       # [has test], [has doc]
    pushevent!,         # [has test], [has doc]
    samplelocalpath     # [has test], [has doc]


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
