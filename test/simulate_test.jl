using PDSampler, Base.Test

# This tests the Simulation object (not the simulate function)

srand(1777)
n = 1000           # n observations
p = 5              # n dimensions (covariates)
X = randn(n,p)+0.1  # feature matrix
w = 5*rand(p)       # true vector of parameters
# observations according to a logistic thresholded to {-1,1}
y = (logistic.(X*w) .> rand(n)) .* 2.0 .- 1.0
# proxy for N*L upper bound
b = sum( mapslices(e->norm(e)^2,X,1) )/4

# DATA MODEL
dm = LogReg(X,y,b)

# GEOMETRY
ns, a    = eye(p), zeros(p)
geom     = Polygonal(ns,a)
nb(x,v)  = nextboundary(geom, x, v)

# GRADIENTS and IPPSampling
gll_cv    = gradloglik_cv(dm, w)
gllstar   = gradloglik(dm, w)
lb        = LinearBound(w, gllstar, dm.b)
nev(x, v) = nextevent_bps(lb, x, v)

# SIMULATION
T           = Inf          # simulation 'time'
maxgradeval = 500000       # max number of gradient evaluations
lref        = 2.           # refreshment rate
x0          = w            # initial point for the trajectory
v0          = randn(dm.p)  # draw velocity from normal distr
v0         /= norm(v0)     # normalise velocity

sim = Simulation(
        x0, v0, T, nev, gll_cv, nb, lref;
        maxgradeval = maxgradeval);

# all the basic ones

xrand = randn(p)
vrand = randn(p)

@test   sim.x0 == x0 &&
        sim.v0 == v0 &&
        sim.T  == T  &&
        (srand(12);sim.nextevent(xrand, vrand).tau)==
        (srand(12);nev(xrand, vrand).tau)
@test   (srand(12);sim.gll(xrand))==(srand(12);gll_cv(xrand)) &&
        sim.nextboundary(xrand, vrand) == nb(xrand, vrand) &&
        sim.lambdaref == lref &&
        sim.algname   == "BPS"
@test   sim.dim         == length(x0) &&
        sim.mass        == eye(0) &&
        sim.blocksize   == 1000 &&
        sim.maxsimtime  == 4e3 &&
        sim.maxsegments == Int(1e6) &&
        sim.maxgradeval == maxgradeval

simf = Simulation(
        x0, v0, T, nev, gll_cv, nb, lref, "bps";
        maxgradeval = maxgradeval);

@test simf.algname == "BPS"

@test_throws AssertionError Simulation(
        x0, v0, T, nev, gll_cv, nb, lref, "bpsasdf";
        maxgradeval = maxgradeval);

sim2 = Simulation(
        x0 = x0, v0 = v0, T = T, nextevent = nev, gradloglik = gll_cv,
        nextboundary = nb, lambdaref = lref; maxgradeval = maxgradeval)

@test   sim2.x0 == x0 &&
        sim2.v0 == v0 &&
        sim2.T  == T  &&
        (srand(12);sim2.nextevent(xrand, vrand).tau)==
        (srand(12);nev(xrand, vrand).tau)
@test   (srand(12);sim2.gll(xrand))==(srand(12);gll_cv(xrand)) &&
        sim2.nextboundary(xrand, vrand) == nb(xrand, vrand) &&
        sim2.lambdaref == lref &&
        sim2.algname   == "BPS"
@test   sim2.dim         == length(x0) &&
        sim2.mass        == eye(0) &&
        sim2.blocksize   == 1000 &&
        sim2.maxsimtime  == 4e3 &&
        sim2.maxsegments == Int(1e6) &&
        sim2.maxgradeval == maxgradeval

@test_throws AssertionError Simulation(
        x0 = x0, v0 = v0, T = T, gradloglik = gll_cv,
        nextboundary = nb, lambdaref = lref; maxgradeval = maxgradeval)
@test_throws AssertionError Simulation(
        x0 = x0, v0 = v0, T = T, nextevent = nev,
        nextboundary = nb, lambdaref = lref; maxgradeval = maxgradeval)
@test_throws AssertionError Simulation(
        x0 = x0, v0 = v0, T = T, nextevent = nev, gradloglik = gll_cv,
        lambdaref = lref; maxgradeval = maxgradeval)
@test_throws AssertionError Simulation(
        x0 = x0, v0 = v0, T = 0.0, nextevent = nev, gradloglik = gll_cv,
        nextboundary = nb, lambdaref = lref; maxgradeval = maxgradeval)
@test_throws AssertionError Simulation(
        x0 = x0, v0 = v0, T = 1.0, nextevent = nev, gradloglik = gll_cv,
        nextboundary = nb, lambdaref = -0.1; maxgradeval = maxgradeval)
@test_throws AssertionError Simulation(
        v0 = v0, T = T, nextevent = nev, gradloglik = gll_cv,
        nextboundary = nb, lambdaref = lref; maxgradeval = maxgradeval)
@test_throws AssertionError Simulation(
        x0 = x0, T = T, nextevent = nev, gradloglik = gll_cv,
        nextboundary = nb, lambdaref = lref; maxgradeval = maxgradeval)
@test_throws AssertionError Simulation(
        x0 = x0, v0 = v0, nextevent = nev, gradloglik = gll_cv,
        nextboundary = nb, lambdaref = lref; maxgradeval = maxgradeval)
