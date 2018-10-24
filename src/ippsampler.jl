export
    NextEvent,
    LinearBound,
    nextevent_bps,
    nextevent_bps_q,
    nextevent_zz

abstract type IPPSamplingMethod end
abstract type Thinning <: IPPSamplingMethod end

"""
    LinearBound

Indicating you have access to a global bound on the eigenvalues of the Hessian.
"""
struct LinearBound <: Thinning
    xstar::Vector{Real}
    gllstar::Vector{Real}  # grad log likelihood around xstar
    b::Real                # N*L where L is Lipschitz constant
    a::Function
    function LinearBound(xstar, gllstar, b)
        new(xstar, gllstar, b,
            (x, v) -> max(dot(-gllstar, v), 0.) + norm(x-xstar) * b)
    end
end

"""
    NextEvent

Object returned when calling a `nextevent_*` type function. The bouncing time
`tau` as well as whether it should be accepted or rejected `dobounce` (set to
true always for analytical sampling). `flipindex` is to support ZZ sampling.
"""
struct NextEvent
    tau::Real # candidate first arrival time
    dobounce::Function # do bounce (if accept reject step)
    flipindex::Int     # flipindex (if ZZ style)
    function NextEvent(tau, dobounce, flipindex)
        new(tau, dobounce, flipindex)
    end
end
NextEvent(tau; dobounce=(g,v)->true, flipindex=-1) =
    NextEvent(tau, dobounce, flipindex)

# -----------------------------------------------------------------------------

"""
    nextevent_bps(g::MvGaussian, x, v)

Return a bouncing time and corresponding intensity in the specific case of a
Gaussian density (for which the simulation of the first arrival time of the
corresponding IPP can be done exactly).
"""
function nextevent_bps(g::MvGaussian,
                       x::Vector{<:Real}, v::Vector{<:Real})
    # precision(g)*x - precision(g)*mu(g) --> precmult, precmu in Gaussian.jl
    a = dot(mvg_precmult(g, x) - mvg_precmu(g), v)
    b = dot(mvg_precmult(g, v), v)
    c = a / b
    e = max(0.0, c) * c # so either 0 or c^2
    tau = -c + sqrt(e + 2randexp() / b)
    return NextEvent(tau)
end

function nextevent_bps(g::PMFGaussian,
                       x::Vector{<:Real}, w::Vector{<:Real})
    # unmasking x, w
    xu, xv = x[1:g.d], x[g.d+1:end]
    wu, wv = w[1:g.d], w[g.d+1:end]
    # precomputing useful dot products
    xuwv, xvwu = dot(xu, wv), dot(xv, wu)
    wuwv, dxw = dot(wu, wv), xuwv + xvwu
    # real root
    t0 = -0.5dxw / wuwv
    # e(x) := ⟨x_u, x_v⟩ - r
    ex = dot(xu, xv) - g.r
    # (quadratic in t) e(x+tw) = e(x)+(⟨wu,xv⟩+⟨wv,xu⟩)t+⟨wu,wv⟩t^2
    p1 = Poly([ex, dxw, wuwv])
    # (linear in t) ⟨xv,wu⟩+⟨xu,wv⟩+2⟨wu,wv⟩t
    p2 = Poly([dxw, 2.0wuwv])
    # ⟨∇E(x+tw),w⟩ is a cubic in t (and the intensity)
    p  = p1 * p2

    rexp = randexp()
    tau = 0.0

    ### CASES (cf. models/pmf.jl)

    Δ = t0^2 - ex / wuwv # discriminant of p1
    if Δ <= 0
        # only single real root (t0)
        # two imaginary roots
        if t0 <= 0
            tau = pmf_caseA(rexp, p)
        else
            tau = pmf_caseB(rexp, p, t0)
        end
    else
        # three distinct real roots,
        tm, tp = t0 .+ sqrt(Δ) * [-1.0,1.0]
        if tp <= 0 # case: |/
            tau = pmf_caseA(rexp, p)
        elseif t0 <= 0 # case: |_/
            tau = pmf_caseB(rexp, p, tp)
        elseif tm <= 0 # case: |\_/
            tau = pmf_caseC(rexp, p, t0, tp)
        else # case: |_/\_/
            tau = pmf_caseD(rexp, p, tm, t0, tp)
        end
    end
    return NextEvent(tau)
end

"""
    nextevent_zz(g::MvGaussian, x, v)

Same as nextevent but for the Zig Zag sampler.
"""
function nextevent_zz(g::MvGaussian,
                      x::Vector{<:Real}, v::Vector{<:Real})
    # # precision(g)*x - precision(g)*mu(g) --> precmult, precmu in Gaussian.jl
    u1 = mvg_precmult(g, x)-mvg_precmu(g)
    u2 = mvg_precmult(g, v)
    taus = zeros(g.p)
    for i in 1:g.p
        ai = u1[i] * v[i]
        bi = u2[i] * v[i]
        ci = ai ./ bi
        ei = max(0.0, ci) * ci
        taus[i] = -ci + sqrt(ei + 2.0randexp() / abs(bi))
    end
    tau, flipindex = findmin(taus)
    return NextEvent(tau, flipindex=flipindex)
end

"""
    nextevent_bps(lb::LinearBound, x, v)

Return a bouncing time and corresponding intensity corresponding to a linear
upperbound described in `lb`.
"""
function nextevent_bps(lb::LinearBound,
                       x::Vector{<:Real}, v::Vector{<:Real})
    a = lb.a(x, v)
    b = lb.b
    @assert a >= 0.0 && b > 0.0 "<ippsampler/nextevent_bps/linearbound>"
    tau = -a / b + sqrt((a / b)^2 + 2randexp() / b)
    lambdabar = a + b * tau
    return NextEvent(tau, dobounce=(g,v)->(rand()<-dot(g, v)/lambdabar))
end

function nextevent_bps_q(gll::Function, x::Vector{<:Real}, v::Vector{<:Real},
                         tref::Real; n=100)

    chi(t) = max(0.0, dot(-gll(x + t * v), v))
    S      = Chebyshev(0.0..tref)
    p      = points(S, n)
    v      = chi.(p)
    f      = Fun(S, ApproxFun.transform(S,v))
    If     = cumsum(f) # integral from 0 to t with t < = tref
    tau    = ApproxFun.roots(If - randexp())
    tau    = (length(tau)>0) ? minimum(tau) : Inf
    return NextEvent(tau)
end
