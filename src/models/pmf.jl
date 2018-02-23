export
    PMFGaussian

struct PMFGaussian
    r::Real   # Rating
    sigma::Real
    d::Int   # Dimension of latent space
end

function gradloglik(g::PMFGaussian, x::Vector{<:Real})
    u = x[1:g.d]
    v = x[g.d+1:end]
    e = dot(u, v) - g.r
    return -e / g.sigma^2 * [v; u]
end

function loglik(g::PMFGaussian, x::Vector{<:Real})
    u = x[1:g.d]
    v = x[g.d+1:end]
    e = dot(u, v) - g.r
    return 0.5(e / g.sigma)^2
end

# --------------------------------------------------------------------------
# HELPER FUNCTIONS FOR THE BPS sampler

function pmf_base(rexp::Real, p::Poly)
    # find solutions to [p(x) = y]
    rs  = roots(p - rexp)
    # we want the first real root > 0.0
    # this requires a bit of subtlety as `roots` is not
    # super stable so there are imaginary numbers trailing
    #
    # start by dumping roots with a negative real parts
    # and roots that are "too imaginary"
    prs = [r for r in rs if (real(r) > 0.0) && (imag(r) / abs(r) < 1e-8) ]
    prs = sort(abs.(prs))
    # now pick the one with the smallest radius
    k   = searchsortedfirst(prs, 0.0)
    return k > length(prs) ? 0.0 : prs[k]
end


function pmf_caseA(rexp::Real, p::Poly)
    # == CASE A ###### intensity:
    #  |   |
    #  |  /
    #  |/
    #  |________________
    ##############################
    # the branch starts on y axis
    # can directly look for root
    ##############################
    return pmf_base(rexp, polyint(p))
end


function pmf_caseB(rexp::T, p::Poly, r::T) where T <: Real
    @assert 0.0 < r
    # == CASE B ###### intensity:
    #   |         |
    #   |        |
    #   |      /
    #   |____r__________
    ##############################################
    # the branch starts after 0 at r
    # move the branch to zero, then add t0 to root
    ##############################################
    return r + pmf_base(rexp, polyint(polyval(p, Poly([r, 1.0]))))
end


function pmf_caseC(rexp::T, p::Poly, t0::T, tp::T) where T <: Real
    @assert 0.0 < t0 < tp
    # == CASE C ###### intensity:
    #   |             |
    #   |\           |
    #   |  \       /
    #   |___t0___tp______
    ##############################################
    # partial hump then flat then branch
    # the integral has a "step", should check on
    # what step we are
    ##############################################
    tau = 0.0
    Ip = polyint(p)
    step = Ip(t0)
    if rexp <= step # solve for 0.0->t0
        tau = pmf_base(rexp, Ip)
    else # solve for after tp
        tau = tp + pmf_base(rexp-step, polyint(polyval(p, Poly([tp, 1.0]))))
    end
    return tau
end


function pmf_caseD(rexp::T, p::Poly, tm::T, t0::T, tp::T
                    ) where T <: Real
    @assert 0.0 <= tm < t0 < tp
    # == CASE D ###### intensity:
    #   |
    #   |       _          |
    #   |     /  \       /
    #   |___tm___t0____tp_____
    ##############################################
    # flat then complete hump then flat then branch
    # the integral has a two "steps", should check
    # on what step we are
    ##############################################
    # start by shifting by tm
    Ip1  = polyint(polyval(p, Poly([tm, 1.0])))
    step = Ip1(t0-tm)
    if rexp <= step # solve between tm and t0
        tau = tm + pmf_base(rexp, Ip1)
    else # solve after tp
        tau = tp + pmf_base(rexp-step, polyint(polyval(p, Poly([tp, 1.0]))))
    end
end
