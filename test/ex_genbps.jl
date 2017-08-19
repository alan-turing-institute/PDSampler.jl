using Base.Test
#@startexample Generalised BPS (Truncated Gaussian)
#=

This example is identical to the Global Bouncy Particle Sampler example with
a truncated gaussian except that it uses the *Generalised* kernel. In
that case, there does not need to be a refreshment step.
See also the corresponding paper: [*Generalized Bouncy Particle Sampler*](https://arxiv.org/pdf/1706.04781.pdf) by Changye Wu and Christian
Robert.

=#
using PDMP
p     = 2
ns, a = eye(p), zeros(p)
geom  = Polygonal(ns, a)

nextbdG(x, v) = nextboundary(geom, x, v)

srand(12)
P1  = randn(p,p)
P1 *= P1'
P1 += norm(P1)/100*eye(p)
C1  = inv(P1); C1 += C1'; C1/=2;
L1  = cholfact(C1)
mu  = zeros(p)+1.
mvg = MvGaussianCanon(mu, P1)

gradllG(x) = gradloglik(mvg, x)

nextevG(x, v) = nextevent_bps(mvg, x, v)
#=
Note the specification of `algname="GBPS"`. Note also the refreshment rate
set to 0.0 (no refreshment).
=#
T    = 1000.0   # length of path generated
lref = 0.0      # rate of refreshment
x0   = mu+L1[:L]*randn(p) # sensible starting point
v0   = randn(p) # starting velocity
v0  /= norm(v0) # put it on the sphere (not necessary)
# Define a simulation
sim = Simulation( x0=x0, v0=v0, T=T, nextevent=nextevG, gradloglik=gradllG,
                  nextboundary=nextbdG, lambdaref=lref, maxgradeval = 10000,algname="GBPS")
#=
The rest is as before:
=#
(path, details) = simulate(sim)

sN = 1000
s  = broadcast(+, mu, L1[:L]*randn(p,sN))
mt = zeros(2)
np = 0
# Sum for all samples in the positive orthan
ss = [s; ones(sN)']
mt = sum(ss[:,i] for i in 1:sN if !any(e->e<0, ss[1:p,i]))
mt = mt[1:p]/mt[end]

#@endexample

# all that's below won't show in doc

@test norm(pathmean(path)-mt)/norm(mt) < 0.05
ess,ns = esspath(path; ns=100)
@test minimum(ess) > 150
