# [Global sampler](@id td-globalsampler)

Link to the source files:

* [`path.jl`](https://github.com/alan-turing-institute/PDSampler.jl/blob/master/src/path.jl)
* [`simulate.jl`](https://github.com/alan-turing-institute/PDSampler.jl/blob/master/src/simulate.jl)

## Path

The `path.jl` file exports a type `Path` and a few functions that apply on such an object.

A `Path` object encapsulates

* `xs` a matrix of size $p \times N_e$ where $p$ is the dimension and $N_e$ are the number of events recorded
* `ts` a vectof of times associated with the events
* `p` the dimension
* `nseg` the number of segments

### Auxiliary type

The type `Segment` is useful to encapsulate the information contained between two events. It contains:

* `ta`, `tb` the times at the two events
* `xa`, `xb` the positions
* `tau` (implicit) the travel time corresponding to the segment or $(t_b-t_a)$
* `v` (implicit) the velocity along the segment

### Auxiliary functions

* `getsegment` retrieves a segment starting at time `Path.ts[j]`
* `samplepath` takes a time or a list of times and returns the position along the `Path` object at those times
* `quadpathpoly` does a simple analytical integration along the path for simple test functions of the form $\varphi(x)=P(x)$ where $P$ is a polynomial (for each dimension)
* `pathmean` computes the mean using `quadpathpoly` where the polynomial is just $x$
* `esspath` computation of the ESS corresponding to a number of samples equally spaced along the path. It tries to achieve a specific ratio and increases the number of samples until it achieves it.

## Simulate

The `simulate.jl` file is the main file for the Global Sampler. It contains one main immutable object which contains all of the parameters of the simulation. This is convenient if multiple parameters need to be tested.

A number of default parameters are pre-encoded but some parameters are essential (starting point, starting velocity, etc).

### Simulate function

This function should mimic rather closely the original paper. Note in the main loop:

```julia
tau = min(bounce.tau, taubd, tauref)
```

corresponds to finding which action needs to be executed (normal bounce, boundary bounce or refreshment). After this one of three branches is executed:

```julia
if tau==bounce.tau
    ...
elseif tau==taubd
    ...
else
    ...
end
```

The first branch corresponds to a bounce against the level set of the log-likelihood, the second to a boundary bounce and third a refreshment.

In the first branch, an explicit call to `bounce.dobounce` checks whether to thin the event or not. If the time is accepted then the velocity is refreshed otherwise the whole loop is ignored.
