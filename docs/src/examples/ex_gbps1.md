# Global BPS (Truncated Gaussian)

(*the code for this example can be found [here](https://github.com/alan-turing-institute/PDMP.jl/blob/master/test/ex_gbps1.jl)*, note that the doc rendered here was automatically generated, if you want to fix it, please do it in the julia code directly*)

In this example we **XXXXXXX**

Start by loading the library:
```julia
using PDMP
```
you will then need to define two elements:
1. a geometry (boundaries)
2. an energy (gradient of the log-likelihood of the target)
At the moment, the package can handle unconstrained geometries and polygonal
domains (**see XXXXX**).
Let's say we want to be constrained to the positive orthan in 2D:
```julia
p = 2
# normal to faces and intercepts
ns, a = eye(p), zeros(p)
geom  = Polygonal(ns, a)
# for a given ray, which boundary does it hit?
nextbd(x, v) = nextboundary(geom, x, v)
```
Here `ns` and `a` are the normals and the intercepts of the faces. The type
`Polygonal` encapsulates the geometry.
The function `nextboundary` returns the next boundary on the current ray
`[x,x+tv]` with `t>0` as well as the time of the hit.

We then need to specify a model: we need to define a function of the form
`gradll(x)` which can return the gradient of the log-likelihood at some point `x`.
Here, let us consider a 2D gaussian.
```julia
# build a valid precision matrix, the cholesky decomposition of
# the covariance matrix will be useful later to build a sensible
# starting point.
srand(12)
P1  = randn(p,p)
P1 *= P1'
P1 += norm(P1)/100*eye(p)
C1  = inv(P1); C1 += C1'; C1/=2;
L1  = cholfact(C1)
mu  = zeros(p)+1.
mvg = MvGaussianCanon(mu, P1)
```
Here, we have defined the gaussian through the "Canonical" representation
(**see XXXXX**) i.e.: by specifying a mean and a precision matrix.

The gradient of the log-likelihood is then given by
```julia
gradll(x) = gradloglik(mvg, x)
```
**Remark**: if you want to implement your own model, you should define your
model in (**XXXXXX**) and make sure it implements a `gradloglik` function.

Next, we need to define the function which can return the first arrival time of
the Inhomogenous Poisson Process (**cf. algorithm**).
Note that you could be using `nextevent_zz` here as well if you wanted to use
the Zig-Zag sampler (and you could implement other kernels as well, see
**HERE XXXXX**).
```julia
nextev(x, v) = nextevent_bps(mvg, x, v)
```
For a Gaussian (and some other simple distributions), this is analytical through
an inversion-like method (**cf. algorithm**).
Another approach is the thinning approach using a bounding intensity.
At the moment thinning with a linear bound is implemented (**cf XXXXX**).

Finally, you need to specify the parameters of the simulation such as the
starting point and velocity, the length of the path generated, the rate of
refreshment and the maximum number of gradient evaluations. (**see discussion**)
```julia
T    = 1000.0   # length of path generated
lref = 2.0      # rate of refreshment
x0   = mu+L1[:L]*randn(p) # sensible starting point
v0   = randn(p) # starting velocity
v0  /= norm(v0) # put it on the sphere (not necessary)
# Define a simulation
sim = Simulation( x0, v0, T, nextev, gradll,
            nextbd, lref ; maxgradeval = 10000)
```
And finally, generate the path and recover some details about the simulation.
```julia
(path, details) = simulate(sim)
```
The `path` object belongs to the type `Path` and can be sampled using
`samplepath`.

A crude test is to check that the estimated mean obtained through quadrature
along the path yields a similar result as a basic Monte Carlo estimator.
```julia
# Building a basic MC estimator
sN = 1000
s  = repmat(mu,1,sN)+L1[:L]*randn(p,sN)
mt = zeros(2)
np = 0
# Sum for all samples in the positive orthan
ss = [s; ones(sN)']
mt = sum(ss[:,i] for i in 1:sN if !any(e->e<0, ss[1:p,i]))
mt = mt[1:p]/mt[end]
```
You can now compare the norm of `mt` to `pathmean(path)` and you will see that
the relative error is below 5%.
