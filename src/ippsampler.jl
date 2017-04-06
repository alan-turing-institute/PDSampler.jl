abstract IPPSamplingMethod
abstract Thinning <: IPPSamplingMethod

# TODO: complete
"""
    LinearBound

Indicating you have access to a global bound on the eigenvalues of the Hessian.
"""
immutable LinearBound <: Thinning
    xstar::Vector{Float}
    gllstar::Vector{Float}  # grad log likelihood around xstar
    b::Float                # N*L where L is Lipschitz constant
    a::Function
    function LinearBound(xstar, gllstar, b)
        new(xstar, gllstar, b,
            (x,v) -> max(dot(-gllstar,v),0.) + norm(x-xstar)*b )
    end
end

# TODO: complete
"""
    NextEvent

Object returned when calling a nextevent_ type function.
"""
immutable NextEvent
    tau::Float         # candidate first arrival time
    dobounce::Function # do bounce (if accept reject step)
    flipindex::Int     # flipindex (if ZZ style)
    function NextEvent(tau; dobounce=(g,v)->true, flipindex::Int=-1)
        new(tau,dobounce,flipindex)
    end
end

# ------------------------------------------------------------------------------

"""
    nextevent_bps(g::MvGaussian, x, v)

Return a bouncing time and corresponding intensity in the specific case of a
Gaussian density (for which the simulation of the first arrival time of the
corresponding IPP can be done exactly).
"""
function nextevent_bps{T<:Vector{Float}}(g::MvGaussian, x::T, v::T)::NextEvent
    # precision(g)*x - precision(g)*mu(g) --> precmult, precmu in Gaussian.jl
    a = dot(mvg_precmult(g,x)-mvg_precmu(g), v)
    b = dot(mvg_precmult(g,v), v)
    c = a/b
    e = max(0.0,c)*c # so either 0 or c^2
    tau = -c + sqrt( e + 2randexp()/b)
    NextEvent(tau)
end

# ### TODO: dev
# function  nextevent_bps{T<:Vector{Float}}(g::PMFGaussian, x::T, w::T)::NextEvent
#     xu,xv = x[1:g.d], x[g.d+1:end]
#     wu,wv = w[1:g.d], w[g.d+1:end]
#     #
#     xuwv = dot(xu,wv)
#     xvwu = dot(xv,wu)
#     wuwv = dot(wu,wv)
#     dxw  = xuwv+xvwu
#     #
#     t0    = -0.5dxw/wuwv
#     ex    = dot(xu,xv)-g.r
#     #
#     # e(x+tw) = e(x)+(<wu,xv>+<wv,xu>)t+<wu,wv>t^2
#     p1 = Poly([ex, dxw, wuwv])
#     # <xv,wu>+<xu,wv>+2<wu,wv>t
#     p2 = Poly([dxw, 2.0wuwv])
#     p  = p1*p2
#     Ip = polyint(p)
#     #
#     delta = t0^2+ex/wuwv
#     rexp = randexp()
#
#     tau = 0.0
#
#     if delta<=0
#         # only single real root (t0)
#         if t0 <= 0
#             tau = pmf_caseA(rexp, Ip)
#         else
#             tau = pmf_caseB(rexp, p, Ip, t0)
#         end
#     else
#         # three distinct real roots, potential bump ==> 0_/\_/ <==
#         tm,tp = t0 + sqrt(t0^2+ex/wuwv)*[-1.0,1.0]
#         if tp <= 0 # ==> 0/ <==
#             # ==CASE A==
#             # no max needed, need to find x s.t. rexp=F(x)
#             intercept = p(0.0)
#             if rexp < intercept
#                 tau = 0.0
#             else
#                 rs  = roots(Ip-rexp)
#                 rrs = [rs[i].re if isapprox(rs[i].im,0.0) for i in 1:length(rs)]
#                 k   = searchsortedfirst(rrs, t0)
#                 tau = rrs[k]
#             end
#         else
#         end
#     end
#
#     ### Cases
#     # Full Hump
#
#
# end

"""
    nextevent_zz(g::MvGaussian, x, v)

Same as nextevent but for the Zig Zag sampler.
"""
function nextevent_zz{T<:Vector{Float}}(g::MvGaussian, x::T, v::T)::NextEvent
    # # precision(g)*x - precision(g)*mu(g) --> precmult, precmu in Gaussian.jl
    u1 = mvg_precmult(g,x)-mvg_precmu(g)
    u2 = mvg_precmult(g,v)
    taus = zeros(g.p)
    for i in 1:g.p
        ai = u1[i]*v[i]
        bi = u2[i]*v[i]
        ci = ai./bi
        ei = max(0.0,ci)*ci
        taus[i] = -ci + sqrt(ei + 2randexp()/abs(bi))
    end
    tau, flipindex = findmin(taus)
    NextEvent(tau,flipindex=flipindex)
end

"""
    nextevent_bps(lb::LinearBound, x, v)

Return a bouncing time and corresponding intensity corresponding to a linear
upperbound described in `lb`.
"""
function nextevent_bps{T<:Vector{Float}}(lb::LinearBound, x::T, v::T)::NextEvent
    a = lb.a(x, v)
    b = lb.b
    @assert a>=0.0 && b>0.0 "<ippsampler/nextevent_bps/linearbound>"
    tau       = -a/b + sqrt((a/b)^2 + 2randexp()/b)
    lambdabar = a + b*tau
    dobounce(g,v) = rand() < -dot(g,v)/lambdabar
    NextEvent(tau, dobounce=dobounce)
end
