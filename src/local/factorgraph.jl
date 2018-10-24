export
    Factor,
    FactorObject,
    FactorGraphStruct,
    FactorGraph,
    assocvariables,
    assocfactors,
    linkedfactors,
    chainstruct,
    chaingraph

"""
    Factor

Encapsulation of a factor, made out of the function that computes the next
bounce, the gradient of the log likelihood associated with the factor and
the index of the factor in the list of factors of a factor graph.
See also: `FactorGraph`, `FactorGraphStruct`.
"""
struct Factor
    nextevent::Function # compute the corresponding first arrival time
    gll::Function       # gradient of log likelihood
    index::Int          # factor index
end

"""
    FactorGraphStruct

Encapsulation of the structure of a factor graph. It is made of a list of
factors each corresponding to a list of variables attached to the factor.
See also: `Factor`, `FactorGraph`.
"""
struct FactorGraphStruct
    flist::Vector{Vector{Int}} # each entry is a list of vars
    # -- implicit
    vlist::Vector{Vector{Int}} # each entry is a list of factors
    nfactors::Int
    nvars::Int
    # -- constructor
    function FactorGraphStruct(flist)
        # finding how many variables there are
        maxvar = 0
        for f in flist
            maxvar = max(maxvar, maximum(f))
        end
        # populating the lists of factors for each variable
        vars = [Vector{Int}(undef, 0) for i in 1:maxvar]
        for (fi, f) in enumerate(flist)
            for e in f
                push!(vars[e], fi)
            end
        end
        new(flist, vars, length(flist), maxvar)
    end
end

"""
    FactorGraph

Encapsulation of a factor graph. It is made out of a structure
(FactorGraphStruct) and an array of factors corresponding to the given
structure.
See also: `Factor`, `FactorGraphStruct`.
"""
struct FactorGraph
    structure::FactorGraphStruct
    factors::Vector{Factor}
end
FactorGraph(fgs::Vector{Vector{Int}}, factors::Vector{Factor}) =
    FactorGraph(FactorGraphStruct(fgs), factors)

# access the variables associated to a factor or vice versa
assocvariables(fgs::FactorGraphStruct, fidx::Int) = fgs.flist[fidx]
assocfactors(fgs::FactorGraphStruct, vidx::Int) = fgs.vlist[vidx]
assocvariables(fg::FactorGraph, fidx::Int) = assocvariables(fg.structure, fidx)
assocfactors(fg::FactorGraph, vidx::Int) = assocfactors(fg.structure, vidx)

"""
    linkedfactors(fgs, fidx)

Return the list of factors that have intersecting variable set with the factor
indexed by `fidx`. It does *not* contain the current factor fidx for efficiency
reason (the current factor can be updated quicker directly from the general
simulation loop).
"""
function linkedfactors(fgs::FactorGraphStruct, fidx::Int)
    lf = Vector{Int}(undef, 0)
    for v ∈ fgs.flist[fidx]
        push!(lf, fgs.vlist[v]...)
    end
    unq = unique(lf)                    # this still contains fidx
    deleteat!(unq, findall((in)(fidx), unq))   # this does not
end
linkedfactors(fg::FactorGraph, fidx::Int) = linkedfactors(fg.structure, fidx)

# ==============================================================================
# Definition of standard FG structures

chainstruct(nvars::Int) = FactorGraphStruct([[i,i+1] for i ∈ 1:nvars-1])
chaingraph(lf::Vector{Factor}) = FactorGraph(chainstruct(length(lf) + 1), lf)

# row based reading (left right then top down)
# gridstruct(nvars::Int) = error("Not implemented yet")
