using PDMP, JLD

srand(1234)

### Define factor graph
# ---------------------

nfac = 3 # number of factors

# description of the likelihood on a factor
mvg             = PDMP.MvGaussianStandard(zeros(2),eye(2))
gll(x)          = PDMP.gradloglik(mvg, x)
nextevent(x, v) = PDMP.nextevent_bps(mvg, x, v)

# all factors have that same likelihood
chainfactor(i) = Factor(nextevent,gll,i)

# assemble into a chain graph
chain = chaingraph([chainfactor(i) for i in 1:nfac])

### Define LocalSimulation
# ------------------------

lambdaref  = .01
maxnevents = 10000
T          = Inf
nvars      = chain.structure.nvars
x0         = randn(nvars)
v0         = randn(nvars)
v0        /= norm(v0)

lsim = LocalSimulation(chain, x0, v0, T, maxnevents, lambdaref)

### Run Local BPS
# ---------------

(all_evlist, details) = simulate(lsim)

### Save results for further inspection
# -------------------------------------

save("res/dex_lbps_gausschain1.jld",
        "evlist", all_evlist.evl,
        "details", details
        )
