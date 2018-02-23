export
    Unconstrained,
    Polygonal,
    nextboundary

abstract type Domain end

"""
    Unconstrained

Signature of domains without boundaries.
"""
struct Unconstrained <: Domain end

"""
    Polygonal

Encapsulate a domain with polygonal boundaries. The boundaries are defined with
the normals to the facets and the intercepts of the corresponding hyperplanes.
"""
struct Polygonal <: Domain
    normals::Matrix{Real} # dim CxD where C=#constraints, D=#features
    intercepts::Vector{Real} # dim C
end

# -----------------------------------------------------------------------------

"""
    nextboundary(ud::Unconstrained, x, v)

Return (NaN, NaN) corresponding to the unconstrained case.
"""
function nextboundary(ud::Unconstrained,
                      x::T, v::T) where T <: Vector{<:Real}
    return (NaN, NaN)
end

"""
    nextboundary(pd::Polygonal, x, v)

Return the time of next boundary hit along the current trajectory when the
domain is polygonal and return the normal of the corresponding boundary.
The current point is `x` and the velocity `v`.
"""
function nextboundary(pd::Polygonal,
                      x::T, v::T) where T <: Vector{<:Real}
    # hitting time along trajectory (x+tv) for a boundary (normal,intercept) is
    # t = intercept/(normal dot v) - (x dot normal)/(normal dot v)
    nsv  = pd.normals * v
    hits = (pd.intercepts - pd.normals * x) ./ nsv
    # hard threshold times with nsv ~ 0 (near-parallel case),
    # remove negative times and times ~ 0 (for numerical stability)
    hits[ map(|, abs.(nsv) .< 1e-10, hits .< 1e-10) ] = Inf
    # get time of hit + index of corresponding boundary
    (t_hit, j) = findmin(hits)
    # return time of hit + normal vector to boundary
    (t_hit, pd.normals[j, :])
end
