export
    NextEvent,
    LinearBound,
    nextevent_bps,
    nextevent_zz

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

# -----------------------------------------------------------------------------

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
        taus[i] = -ci + sqrt(ei + 2.0randexp()/abs(bi))
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
