using PDSampler
using Random

# ------------------------------------------------------------------------------
# Quick version of functions (for testing) (drawn from some of SJV's old code)
function q_nextboundary(ns, a, x, v)
    nconstr = length(a)
    bts = a./(ns*v) - ns*x./(ns*v)
    mask1 = Bool[isapprox((ns*v)[i], 0.0, atol=1.e-10)  for i=1:nconstr]
    mask2 = Bool[isapprox(bts[i], 0.0, atol=1.e-10)  for i=1:nconstr]
    for (i, e) ∈ enumerate(mask1)
        e && (bts[i] = Inf)
    end
    for (i, e) ∈ enumerate(mask2)
        e && (bts[i] = Inf)
    end
    bts+=(bts.<0.0)*Inf
    (t,j)=findmin(bts)
    (t,vec(ns[j,:]))
end
# ------------------------------------------------------------------------------

Random.seed!(123)

p = 50
c = 10
normals = randn(c,p)
intercepts = randn(c)

unconstr = Unconstrained()
polyconstr = Polygonal(normals,intercepts)

x = randn(p)
v = randn(p)

res = nextboundary(unconstr, x, v)

@test reduce(*, map(isnan, res))

(t_b, n_b) = q_nextboundary(normals, intercepts, x, v)
(t, n)     = nextboundary(polyconstr, x, v)

@test abs(t_b - t) <= 1e-12 && norm(n_b - n) <= 1e-12
