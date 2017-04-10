# PDMP

[![Build Status](https://travis-ci.org/alan-turing-institute/PDMP.jl.svg?branch=master)](https://travis-ci.org/alan-turing-institute/PDMP.jl)

[![codecov.io](http://codecov.io/github/alan-turing-institute/PDMP.jl/coverage.svg?branch=master)](http://codecov.io/github/alan-turing-institute/PDMP.jl?branch=master)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://alan-turing-institute.github.io/PDMP.jl/latest)

The aim of this package is to provide a flexible and expandable framework for the methods discussed in

* Alexandre Bouchard-Côté, Sebastian J. Vollmer and Arnaud Doucet, [The Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain Monte Carlo Method](https://arxiv.org/abs/1510.02451) (preprint, 2015)
* Joris Bierkens, Alexandre Bouchard-Côté, Arnaud Doucet, Andrew B. Duncan, Paul Fearnhead, Gareth Roberts and Sebastian J. Vollmer, [Piecewise Deterministic Markov Processes for Scalable Monte Carlo on Restricted Domains](https://arxiv.org/pdf/1701.04244.pdf) (preprint, 2017)
* Joris Bierkens, Paul Fearnhead and Gareth Roberts, [The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data](https://arxiv.org/pdf/1607.03188.pdf) (preprint, 2016)

## Package under construction

This package is currently under construction.
If you have questions / requests, please send an email to `tlienart σ turing.ac.uk` or open an issue.
If you want to try it out, you should be able to add it using the git url:

#### Install from Github for the first time
```julia
Pkg.clone("git://github.com/alan-turing-institute/PDMP.jl.git")
```
#### Run tests
```julia
Pkg.test("PDMP")
```
#### Update to latest Github version
```julia
Pkg.update("PDMP")
```

**Note about code coverage**: two significant files (`src/simulate.jl` and `src/local/simulate.jl`) are hard to "cleanly" check in the usual CS sense. The guarantees are all in terms of infinite computational time which obviously we do not have.
Some of the code will in the future be pulled out of the main function to have controllable chunks.

## Examples

The folder `examples` contains a number of curated examples for either the global PDMP sampler or the local BPS.
Examples such as `ex_bps_mvg1.jl` use the global BPS approach, `mvg1` refers to the case (here a multivariate gaussian).
Examples such as `ex_lbps_gausschain1` use the local BPS approach, `gausschain1` refers to a Markov chain with Gaussian potentials.
See for example `ex_lbps_gausschain1` which runs the local BPS sampler on a simple chain with standard gaussian factors.

You can have a look at the pre-run [notebook](https://github.com/alan-turing-institute/PDMP/blob/master/examples/analysis_global.ipynb). (Global BPS)
You can do the same with [this notebook](https://github.com/alan-turing-institute/PDMP/blob/master/examples/analysis_local.ipynb). (Local BPS)

The examples produce `.jld` files which are/can then be analysed in one of the `analysis` notebooks (or your adapted version theoref).

## How to use this framework

### Global BPS

(We follow here step by step the example `ex_bps_mvg3.jl`)

Start by loading the library:

```julia
using PDMP
```

you will then need to define two things:

1. a geometry (boundaries)
2. an energy (gradient of log-likelihood)

At the moment, the package can handle unconstrained geometries and polygonal domains.
Let's say we want to be constrained to the positive orthan in 2D:

```julia
p                 = 2
ns, a             = eye(p), zeros(p)
geom              = PDMP.Polygonal(ns,a)
nextboundary(x,v) = PDMP.nextboundary(geom, x, v)
```

Here `ns` and `a` are the normals and the intercepts of the facets.
The type `Polygonal` encapsulates that geometry.
The function `nextboundary` returns the next boundary on the current ray `[x,x+tv]` with `t>0` as well as the time when that hit will happen.

We then need to specify a model, in particular we need to define a function of the form `gll(x)` which can return the gradient of the log-likelihood at `x`.
Here let us consider a 2D gaussian for simplicity.

```julia
P1  = randn(p,p)
P1 *= P1'
mu  = zeros(p)+1.
mvg = PDMP.MvGaussianCanon(mu, P1)
```

Here we have defined the gaussian through the "Canonical" representation (see `src/models/gaussian.jl` for others) i.e.: by specifying a mean and a precision matrix.

The gradient of the log-likelihood is then given by

```julia
gll(x) = PDMP.gradloglik(mvg, x)
```

**Remark**: if you want to implement your own model, you should define your model in `src/models/yourmodel.jl` and make sure it implements a `gradloglik` function.

Next we need to define the function which can return the first arrival time of the IPP (see algorithm). Note that you could be using `nextevent_zz` here as well if you wanted the Zig-Zag sampler. See `src/ippsampler.jl` for more.

```julia
nextevent(x, v) = PDMP.nextevent_bps(mvg, x, v)
```

For a Gaussian (and some other simple functions), this is analytical through an inversion-like method (see BPS paper).
Another approach is the thinning approach using a bounding intensity.
At the moment thinning with a linear bound is implemented (cf. `src/ippsampler.jl`).

Finally, you need to specify the parameters of the simulation such as the starting point, velocity, length of the path generated, rate of refreshment and maximum number of gradient evaluations.

```julia
T    = 1000.0   # length of path generated
lref = 2.0      # rate of refreshment
x0   = randn(p) # starting point
v0   = randn(p) # starting velocity
v0  /= norm(v0) # put it on the sphere (not necessary)

sim = PDMP.Simulation( x0, v0, T, nextevent, gll,
            nextboundary, lref ; maxgradeval = 10000)
```

And finally, generate the path and recover some details about the simulation.

```julia
(path, details) = PDMP.simulate(sim)
```

The `path` object belongs to the type `Path` and can be sampled using `samplepath` (see `src/path.jl`).

### Local BPS

(We follow here step by step the example `ex_lbps_gausschain1.jl`)

The approach to using the local BPS is much the same as for the global one except that you need to specify a `FactorGraph`. That object will contain the structure of the factor graph (which factor is connected to which variables) as well as the list of all factors (which have a `gll` and `nextevent` since each factor can be seen individually as a small BPS).

Let's declare a chain of bivariate gaussians:

```julia
nfac = 3 # number of factors

mvg             = PDMP.MvGaussianStandard(zeros(2),eye(2))
gll(x)          = PDMP.gradloglik(mvg, x)
nextevent(x, v) = PDMP.nextevent_bps(mvg, x, v)

# all factors have that same likelihood
chainfactor(i) = Factor(nextevent,gll,i)

# assemble into a chain graph
chain = chaingraph([chainfactor(i) for i in 1:nfac])
```

This is a simple graph with a known structure so that it's already defined through the `chaingraph` function (in `src/local/factorgraph.jl`)
For an arbitrary graph, you would need to provide two things:

1. the structure of the factor graph: a list of list where each element corresponds to a factor and the corresponding list contains the indices of the variables attached to that factor
2. the list of factors

The rest is very similar to the global BPS:

```julia
lambdaref  = .01
maxnevents = 10000
T          = Inf
nvars      = chain.structure.nvars
x0         = randn(nvars)
v0         = randn(nvars)
v0        /= norm(v0)

lsim = LocalSimulation(chain, x0, v0, T, maxnevents, lambdaref)

(all_evlist, details) = simulate(lsim)
```

The `all_evlist` object contains a list of `EventList` corresponding to the what happened on each of the factors.
It can also be sampled using `samplelocalpath` (cf. `src/local/event.jl`).

# Technical doc

**Comment**: this is bound to become a stand-alone document.

## Understanding the code

### Type aliases and Union types

* `Int` for `Int64` -- in `PDMP.jl`
* `Float` for `Float64` -- in `PDMP.jl`
* `AllowedVarType` for `Union{Float, Vector{Float}}` (variable types in the local BPS) -- in `local/event.jl`
* `AllowedTimeType` for `Union{Vector{Float}, LinSpace{Float}}` (format of collections of times that can be passed to extract samples from events) -- in `path.jl`

### Abstract types

* `MvGaussian` for multivariate gaussians -- in `models/mvgaussian.jl`
* `Domain` for domain description (in global BPS) -- in `geometry.jl`
* `IPPSamplingMethod` for sampling methods of an IPP -- in `ippsampler.jl`
* `Thinning <: IPPSamplingMethod` for thinning methods of an IPP -- in `ippsampler.jl`

### Specific types (global)
**geometry**
* `Unconstrained <: Domain` (**immutable**) defines a domain without boundary (signature only)
* `Polygonal <: Domain` (**immutable**) defines a domain with affine boundaries determined by `normals` and `intercepts`

**ippsampler**
* `LinearBound <: Thinning` (**immutable**) thinning method using a linear bound following a uniform bound on the eigenvalues of the Hessian
* `NextEvent` object returned when sampling from the IPP (contains a bouncing time, whether to bounce or not and a flipindex (ZZ case))

**path**
* `Path` container for a path of the global sampling (stores list of corners and times)
* `Segment` (**immutable**) container to store the two corners of a linear segment, makes it easier to sample from a `Path`

**simulate**
* `Simulation` container for the parameters of a global simulation

### Model types (global)

* `LogReg` (**immutable**) model for a logistic regression

### Specific types (local)

**event**
* `Event{T<:AllowedVarType}` (**immutable**) a triple `(x,v,t)` corresponding to an event for one of the node.
* `EventList{T<:AllowedVarType}` list of events (stored as lists of `xs`, `vs`, `ts` in order to be traversed more effectively than a `Vector{Event}`).
* `AllEventList` container for the `EventList` of all the nodes

**factorgraph**
* `Factor` (**immutable**) attaches the `nextevent` function (sampling from IPP), the gradient of the corresponding log-likelihood (`gll`) and an index
* `FactorGraphStruct` (**immutable**) contains the pattern of connections between nodes via `flist` and `vlist` (storing what variable is attached to what factor and vice-versa) `nfactors` and `nvars` keep track of the number of factors and the number of nodes (variables)
* `FactorGraph` (**immutable**) container for all the `Factor` and the `FactorGraphStruct`

**simulate**
* `LocalSimulation` (**immutable**) container for the parameters of a local simulation (the factor graph, initial positions, etc)
