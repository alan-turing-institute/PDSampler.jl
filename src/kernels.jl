export
    reflect_bps!,
    reflect_zz!,
    refresh_global!,
    refresh_restricted!,
    refresh_partial!,
    reflect_gbps

# ------------------------------------------------------------------------------
# BPS KERNELS

"""
    reflect_bps!(n, v)

BPS specular reflection (in place) of a velocity `v` against a plane defined by
its normal `n`.
"""
function reflect_bps!(n::Vector{<:Real}, v::Vector{<:Real})
    v .-= 2.0dot(n, v) * n / dot(n, n)
    return v
end

"""
    reflect_gbps(n, v)

(Generalized) BPS specular reflection - flip normal component and resample
orthogonal https://arxiv.org/abs/1706.04781
"""
function reflect_gbps(n::Vector{<:Real}, v::Vector{<:Real})
    v2 = randn(length(n))
    v2 .-= dot(n, v2 + v) * n / dot(n, n)
    return v2
end

# """
#     reflect_gbps_sphere!(n, v)
#
# (Generalized) BPS specular reflection flip normal component and resample orthogonal https://arxiv.org/abs/1706.04781
#
# """
# function reflect_gbps_sphere!{T<:Vector{Real}}(n::T, v::T)::T
#     p = dot(n, v)*n/norm(n)^2
#     v -= dot(n, v)*n/norm(n)^2
#     v2 = rnorm(length(d))
#     v2 -= dot(n, v2)*n/norm(n)^2
#     v=-p+v2*norm(v)/norm(v2)
# end

"""
    reflect_bps!(n, v, mass)

BPS specular reflection (in place) of a velocity `v` against a plane defined by
its normal `n` and a mass matrix `mass`.
"""
function reflect_bps!(n::Vector{<:Real}, v::Vector{<:Real},
                      mass::Matrix{<:Real})
    v .-= 2.0dot(n, v) * (mass * n) / dot(mass * n, n)
    return v
end

# ----------------------------------------------------------------------------
# ZZ KERNELS

"""
    reflect_zz!(mask, v)

ZigZag reflection (in place) of v according to a mask on its indeces `mask`.
"""
reflect_zz!(mask::Vector{Int}, v::Vector{<:Real}) = (v[mask] .*= -1.0; v)

reflect_zz!(mask::Int, v::Vector{<:Real}) = reflect_zz!([mask], v)

# --------------------------------------------------------------------------
# Refreshment kernels

refresh_global!(v::Vector{<:Real}) = (v .= randn(length(v)); v)

function refresh_restricted!(v::Vector{<:Real})
    v  = refresh_global!(v)
    v /= norm(v)
    return v
end

function refresh_partial!(v::Vector{<:Real}, beta::Beta{<:Real})
    # sample a vector with norm 1
    w = randn(length(v))
    w /= norm(w)
    # sample an angle
    angle = rand(beta) * 2 * pi
    # scale it by tan(angle):
    #   /|
    # /__|
    # square triangle with base 1,
    # (angle) on left => opposite side=tan(angle)
    w *= tan(angle)
    # take vector on the cone
    v += w
    # normalise
    v /= norm(v)
end
refresh_partial!(beta::Beta{<:Real}) = v -> refresh_partial!(v, beta)
