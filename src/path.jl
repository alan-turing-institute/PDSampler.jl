export
    Path,
    samplepath,
    quadpathpoly,
    pathmean,
    esspath

"""
    AllowedTimeType

Syntactic sugar for union type of Vector{Float} and LinSpace{Float} (types
accepted for the `samplelocalpath` function).
"""
AllowedTimeType = Union{Vector{Float}, Range{Float}}


"""
    Path

Type to store a path simulated via PDMP sampling. It stores the corners and the
times.
"""
type Path
    xs::Matrix{Float}  # corners (dim: nfeatures * nsteps)
    ts::Vector{Float}  # times (dim: nsteps)
    # -- implicit
    p::Int      # nfeatures
    nseg::Int   # number of segments
    Path(xs,ts) = new(xs,ts,size(xs,1), length(ts)-1)
end

"""
    Segment

Type to store the information relative to a single segment within a path in the
global BPS.
"""
immutable Segment
    ta::Float
    tb::Float
    xa::Vector{Float}
    xb::Vector{Float}
    # -- implicit
    tau::Float
    v::Vector{Float}
    function Segment(ta,tb,xa,xb)
        @assert tb>ta "times of subsequent events should be distinct"
        new(ta,tb,xa,xb,tb-ta,(xb-xa)/(tb-ta))
    end
end

"""
    getsegment(p, j)

Retrieves segment j where j=1 corresponds to the first segment from the initial
time to the first corner.
"""
getsegment(p::Path, j::Int) = Segment(p.ts[j],p.ts[j+1],p.xs[:,j],p.xs[:,j+1])

"""
    samplepath(p, t)

Sample the piecewise linear trajectory defined by the corners `xs` and the
times `ts` at given time `t`. The matrix returned has dimensions `p*length(t)`.
"""
function samplepath(p::Path, t::AllowedTimeType)::Matrix{Float}
    @assert t[1] >= p.ts[1] && t[end] <= p.ts[end]
    samples = zeros(p.p,length(t))
    #
    ti  = t[1]
    j   = searchsortedlast(p.ts,ti)
    seg = getsegment(p,j)
    i   = 1
    while i <= length(t)
        ti = t[i]
        if ti <= seg.tb
            samples[:,i] = seg.xa + (ti-seg.ta) * seg.v
            # increment
            i += 1
        else
            j  += searchsortedlast(p.ts[j+1:end],ti)
            # safeguard for out of bound: point to last segment
            j   = (j+1>=length(p.ts)) ? length(p.ts)-1 : j
            seg = getsegment(p,j)
        end
    end
    samples
end
samplepath(p::Path, t::Float) = samplepath(p,[t])[:,1]

"""
    quadpathpoly(path, poly)

Integrate a polynomial on the elements along a path. So if we are sampling in Rd
and have phi(X)=Poly(xk), then we can integrate the polynomial along the path.
Note that we can't do the same thing if phi(X) mixes 2 or more components.
"""
function quadpathpoly(path::Path, pol::Poly)::Vector{Float}
    dim  = length(path.xs[:,1])
    res  = zeros(dim)
    nseg = length(path.ts)-1
    T    = path.ts[end]
    for segidx = 1:nseg-1
        # for each segment in the path
        segment = getsegment(path, segidx)
        # get the characteristics
        tau = segment.tau
        xa  = segment.xa
        v   = segment.v
        # integrate along that segment
        for d = 1:dim
            res[d] += polyint(polyval(pol,Poly([xa[d],v[d]])))(tau)
        end
    end
    # last segment
    lastsegment = getsegment(path, nseg)
    # characteristics
    tau = T - lastsegment.ta
    xa  = lastsegment.xa
    v   = lastsegment.v
    # integrate along that segment
    for d = 1:dim
        res[d] += polyint(polyval(pol,Poly([xa[d],v[d]])))(tau)
    end
    res/T
end
pathmean(path::Path) = quadpathpoly(path, Poly([0.0,1.0]))

"""
    esspath(path, ns, rtol, limitfrac)

Increase the number of sample points on the path. The hope is that the ESS for
each dimension divided by the number of samples becomes larger than `limfrac`.
"""
function esspath(path::Path; ns::Int=1000, rtol::Float=0.1, limfrac::Float=0.2,
                maxns::Int=100000)::Tuple{Vector{Float},Int}

    Tp = 0.999 * path.ts[end]
    gg = linspace(0, Tp, ns)

    spath  = samplepath(path, gg)
    oldess = [ess(spath[i,:]) for i in 1:path.p]
    flag   = any( oldess./ns .> limfrac )
    curess = flag ? similar(oldess) : oldess

    while flag && (ns < maxns/2)
        ns    *= 2
        gg     = linspace(0, Tp, ns)
        spath  = samplepath(path, gg)
        curess = [ess(spath[i,:]) for i in 1:path.p]

        # cond1 : check whether the ESS/NS is big enough
        # cond2 : check whether the ESS changed significantly
        flag   = any( curess./ns .> limfrac) ||
                 any( (abs.(curess - oldess) ./ oldess) .> rtol )
        oldess = curess
    end
    (curess, ns)
end
