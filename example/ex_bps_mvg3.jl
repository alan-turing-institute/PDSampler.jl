using PDMP, JLD


### Compute mean of truncated gaussian
# ------------------------------------

srand(234)

p   = 2
P1  = randn(p,p)
P1 *= P1'
mu  = zeros(p)+1.
C1  = inv(P1)
C1 += C1'
C1 /= 2

sN  = 10000
spl = rand(Distributions.MvNormal(C1),sN)
for i in 1:p
    spl[i,:] += mu[i]
end

splp_m = zeros(2)
np     = 0
for i in 1:sN
    if !any(e->e<0, spl[:,i])
        splp_m += spl[:,i]
        np     += 1
    end
end
trunc_mean = splp_m / np

### Define a MVG
# --------------

T    = 1000.0
lref = 2.

srand(234)

p   = 2
P1  = randn(p,p)
P1 *= P1'
mu  = zeros(p)+1.

# DATA MODEL
mvg = PDMP.MvGaussianCanon(mu, P1)

# GEOMETRY
ns, a             = eye(p), zeros(p)
geom              = PDMP.Polygonal(ns,a)
nextboundary(x,v) = PDMP.nextboundary(geom, x, v)

# GRADIENTS + IPP sampler
gll(x)          = PDMP.gradloglik(mvg, x)
nextevent(x, v) = PDMP.nextevent_bps(mvg, x, v)

# Simulation
x0  = mu            # initial point
v0  = randn(mvg.p)  # draw velocity
v0 /= norm(v0)      # normalise velocity

sim = PDMP.Simulation(x0, v0, T, nextevent, gll,nextboundary, lref;
                        maxgradeval=10000)

### Run Global BPS
# ----------------
(path, details) = PDMP.simulate(sim)

### Save results for further inspection
# -------------------------------------

save("res/dex_bps_mvg3.jld",
        "path", path,
        "details", details,
        "mvg", mvg,
        "trunc_mean", trunc_mean)
