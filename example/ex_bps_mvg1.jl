using PDMP, JLD

srand(1235)

### Define a MVG
# --------------

p   = 10
P1  = randn(p,p)
P1 *= P1'
mu  = zeros(p)

# DATA MODEL
mvg = PDMP.MvGaussianCanon(mu, P1)

# GEOMETRY
geom = PDMP.Unconstrained()
nextboundary(x, v) = PDMP.nextboundary(geom, x, v)

# GRADIENTS and IPPSampling
gll(x) = PDMP.gradloglik(mvg, x)
nextevent(x, v) = PDMP.nextevent_bps(mvg, x, v)

# SIMULATION

T, lref     = Inf, 2.
maxgradeval = 10000

x0  = mu            # initial point
v0  = randn(mvg.p)  # draw velocity
v0 /= norm(v0)      # normalise velocity

sim = PDMP.Simulation(x0, v0, T, nextevent, gll, nextboundary, lref;
                        maxgradeval=maxgradeval)

### Run Global BPS
# ----------------
(path, details) = PDMP.simulate(sim)

### Save results for further inspection
# -------------------------------------

save("res/dex_bps_mvg1.jld",
        "path", path,
        "details", details,
        "mvg", mvg
        )
