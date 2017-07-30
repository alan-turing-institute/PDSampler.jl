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
function reflect_bps!{T<:Vector{Float}}(n::T, v::T)::T
    v -= 2.0dot(n, v)*n/norm(n)^2
    v
end

"""
    reflect_gbps(n, v)

(Generalized) BPS specular reflection - flip normal component and resample orthogonal https://arxiv.org/abs/1706.04781
"""
function reflect_gbps{T<:Vector{Float}}(n::T, v::T)::T
    v2  = randn(length(n))
    v2 -= dot(n, v2 + v)*n/norm(n)^2
    v2
end

# """
#     reflect_gbps_sphere!(n, v)
#
# (Generalized) BPS specular reflection flip normal component and resample orthogonal https://arxiv.org/abs/1706.04781
#
# """
# function reflect_gbps_sphere!{T<:Vector{Float}}(n::T, v::T)::T
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
function reflect_bps!{T<:Vector{Float}}(n::T, v::T, mass::Matrix{Float})::T
    v -= 2.0dot(n, v)*(mass*n)/dot(mass*n,n)
    v
end

# ------------------------------------------------------------------------------
# ZZ KERNELS

"""
    reflect_zz!(mask, v)

ZigZag reflection (in place) of v according to a mask on its indeces `mask`.
"""
function reflect_zz!{T<:Vector{Float}}(mask::Vector{Int}, v::T)::T
    v[mask] *= -1.0
    v
end
reflect_zz!{T<:Vector{Float}}(mask::Int, v::T)::T = reflect_zz!([mask], v)

# ------------------------------------------------------------------------------
# Refreshment kernels

function refresh_global!{T<:Vector{Float}}(v::T)::T
    v = randn(length(v))
    v
end
function refresh_restricted!{T<:Vector{Float}}(v::T)::T
    v  = refresh_global!(v)
    v /= norm(v)
end
function refresh_partial!{T<:Vector{Float}}(v::T, beta::Beta{Float})::T
    # sample an angle
    angle = rand(beta) * 2 * pi
    # sample a vector with norm 1
    w  = randn(length(v))
    w /= norm(w)
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
refresh_partial!(beta::Beta{Float}) = v->refresh_partial!(v,beta)
