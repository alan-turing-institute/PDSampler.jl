using PDMP, Base.Test, Polynomials

d = 3
u = randn(d)
v = randn(d)
r = 3
s = 0.1

pmfg = PMFGaussian(r, s, d)

e  = (dot(u,v) - r)
lb = 0.5(e/s)^2
l  = loglik(pmfg, [u;v])

@test abs(l-lb) <= 1e-10

gub = - e * v / s^2
gvb = - e * u / s^2
g   = gradloglik(pmfg, [u;v])

@test   norm(gub-g[1:d]) <= 1e-10 &&
        norm(gvb-g[d+1:end]) <= 1e-10

### test pmf_base corner cases
rex = 0.1
r   = -0.5+sqrt(0.75)im
p   = poly([1.2+1e-9im,-1.7,r,conj(r)]) + rex
@test isapprox(PDMP.pmf_base(rex, p), 1.2)

# the tests have crude accuracy (because the benchmark is drawn from visual
# inspections with Grapher). This doesn't matter much.

# CASE A # p = (s+1)(s^2+1)
p = Poly([1.0,1.0])*Poly([1.0,0.0,1.0])

@test abs(PDMP.pmf_caseA(2.0, p) - 0.9782) <= 1e-3

# CASE B # p = (s-1)(s^2+1)
p = Poly([-1.0,1.0])*Poly([1.0,0.0,1.0])

@test abs(PDMP.pmf_caseB(2.0, p, 1.0) - 2.0165) <= 1e-3

# CASE C # p (s+1)(s-1)(s-2)
p = Poly([1.0,1.0])*Poly([-1.0,1.0])*Poly([-2.0,1.0])

@test abs(PDMP.pmf_caseC(0.5014, p, 1.0, 2.0) - 0.2761) <= 1e-3
@test abs(PDMP.pmf_caseC(1.4090, p, 1.0, 2.0) - 2.3968) <= 1e-3

# CASE D # p = (s-1)(s-2)(s-3)
p = Poly([-1.0,1.0])*Poly([-2.0,1.0])*Poly([-3.0,1.0])

@test abs(PDMP.pmf_caseD(0.1398, p, 1.0, 2.0, 3.0) - 1.4978) <= 1e-3
@test abs(PDMP.pmf_caseD(0.8405, p, 1.0, 2.0, 3.0) - 3.5928) <= 1e-3
