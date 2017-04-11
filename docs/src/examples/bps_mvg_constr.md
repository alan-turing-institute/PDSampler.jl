# Truncated Multivariate Gaussian

(*we follow here step by step [this example](https://github.com/alan-turing-institute/PDMP.jl/blob/master/example/bps_mvg_constr.jl)*)

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
