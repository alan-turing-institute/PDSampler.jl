immutable PMFGaussian
    r::Float       # Rating
    sigma::Float
    d::Int   # Dimension of latent space
end

function gradloglik(g::PMFGaussian, x::Vector{Float})::Vector{Float}
    u = x[1:g.d]
    v = x[g.d+1:end]
    e = dot(u,v)-g.r
    -e/g.sigma^2 * [v; u]
end

function loglik(g::PMFGaussian, x::Vector{Float})::Float
    u = x[1:g.d]
    v = x[g.d+1:end]
    e = dot(u,v) - g.r
    0.5(e/g.sigma)^2
end

# --------------------------------------------------------------------------
# HELPER FUNCTIONS FOR THE BPS sampler

function pmf_base(rexp::Float, poly::Poly)::Float
    # find solutions to [Poly(x) = y]
    rs  = roots(poly-rexp)
    # pick the one that works for us (first real one after 0.0)
    rrs = [real(rs[i]) for i in 1:length(rs) if isapprox(imag(rs[i]),0.0)]
    rrs = sort(rrs)
    k   = searchsortedfirst(rrs, 0.0)
    rrs[k]
end
function pmf_caseA(rexp::Float, p::Poly)::Float
    # == CASE A ###### intensity:
    #  |   |
    #  |  /
    #  |/
    #  |________________
    ##############################
    # the branch starts on y axis
    # can directly look for root
    ##############################
    pmf_base(rexp, polyint(p))
end
function pmf_caseB{T<:Float}(rexp::T, p::Poly, r::T)::T
    # == CASE B ###### intensity:
    #   |         |
    #   |        |
    #   |      /
    #   |____/__________
    ##############################################
    # the branch starts after 0 at t0
    # move the branch to zero, then add t0 to root
    ##############################################
    r + pmf_base(rexp, polyint(polyval(p, Poly([r,1.0]))))
end
function pmf_caseC{T<:Float}(rexp::T, p::Poly, t0::T, tp::T)::T
    # == CASE C ###### intensity:
    #   |             |
    #   |\           |
    #   |  \       /
    #   |___\____/______
    ##############################################
    # partial hump then flat then branch
    # the integral has a "step", should check on
    # what step we are
    ##############################################
    tau  = 0.0
    Ip   = polyint(p)
    step = Ip(t0)
    if rexp <= step # solve for 0.0->t0
        tau = pmf_base(rexp, Ip)
    else # solve for after tp
        tau = tp + pmf_base(rexp-step, polyint(polyval(p, Poly([tp,1.0]))))
    end
    tau
end
function pmf_caseD{T<:Float}(rexp::T, p::Poly, tm::T, t0::T, tp::T)::T
    # == CASE D ###### intensity:
    #   |
    #   |       _          |
    #   |     /  \       /
    #   |___/_____\____/_____
    ##############################################
    # flat then complete hump then flat then branch
    # the integral has a two "steps", should check
    # on what step we are
    ##############################################
    # start by shifting by tm
    Ip1  = polyint(polyval(p, Poly([tm,1.0])))
    step = Ip1(t0-tm)
    if rexp <= step # solve between tm and t0
        tau = tm + pmf_base(rexp, Ip1)
    else # solve after tp
        tau = tp + pmf_base(rexp-step, polyint(polyval(p, Poly([tp,1.0]))))
    end
end
