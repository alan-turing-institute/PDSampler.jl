using PDMP, JLD

srand(1777)
n = 20000           # n observations
p = 40              # n dimensions (covariates)
X = randn(n,p)+0.1  # feature matrix
w = 5*rand(p)       # true vector of parameters
# observations according to a logistic thresholded to {-1,1}
y = (PDMP.logistic.(X*w) .> rand(n)) .* 2.0 .- 1.0
# proxy for N*L upper bound
b = sum( mapslices(_->norm(_)^2,X,1) )/4

# DATA MODEL
dm = PDMP.LogReg(X,y,b)

# GEOMETRY
ns, a             = eye(p), zeros(p)
geom              = PDMP.Polygonal(ns,a)
nextboundary(x,v) = PDMP.nextboundary(geom, x, v)

# GRADIENTS and IPPSampling
gll_cv          = PDMP.gradloglik_cv(dm, w)
gllstar         = PDMP.gradloglik(dm, w)
lb              = PDMP.LinearBound(w, gllstar, dm.b)
nextevent(x, v) = PDMP.nextevent_bps(lb, x, v)

# SIMULATION
T           = Inf          # simulation 'time'
maxgradeval = 500000       # max number of gradient evaluations
lref        = 2.           # refreshment rate
x0          = w            # initial point for the trajectory
v0          = randn(dm.p)  # draw velocity from normal distr
v0         /= norm(v0)     # normalise velocity

sim = PDMP.Simulation(x0, v0, T, nextevent, gll_cv, nextboundary, lref; maxgradeval = maxgradeval);

(path, details) = PDMP.simulate(sim);


save("res/dex_bps_logreg1.jld",
        "path", path,
        "details", details,
        "dm", dm,
        "w", w)
