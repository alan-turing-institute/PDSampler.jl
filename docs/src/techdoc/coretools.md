# [Core tools](@id td-coretools)

Link to the source files:

* [`geometry.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/geometry.jl)
* [`ippsampler.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/ippsampler.jl)
* [`kernels.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/kernels.jl)

## Geometry

The `geometry.jl` code exports the following **immutable** types:

* `Unconstrained <: Domain `,
* `Polygonal <: Domain`.

and a function

* `nextboundary`

Any geometry type needs to have a `nextboundary` function associated with it of the form `nextboundary(g, x, v)` where `g` refers to a geometry object, `x` to the position and `v` to the velocity. The aim of that function is to:

* return a time of first hit $\tau_h$ (following the ray $x+vt$, for $t>0$),
* return the normal to the boundary at that hitting point.

So for example, the `Unconstrained` is an empty type and the associated `nextboundary` function just returns `(NaN, NaN)` indeed a boundary will never be crossed and there is no normal. These `NaN` are processed by calling functions.

The `Polygonal` domain requires the definition of the normals and the intercepts. The `nextboundary` function is pretty simple (intersection of lines). Note a few tricks for numerical stability:

* near parallel case can be ignored (the crossing time will be huge compared to bouncing time or refreshment time and therefore we can just ignore it)
* negative times can be ignored (we're only going forward)
* times that are very close to zero are ignored (means that we are currently already very close to the boundary meaning that we will bounce away)


## Sampling from an IPP

The `ippsampler.jl` exports the following **immutable** types:

* `NextEvent`,
* `LinearBound <: Thinning <: IPPSamplingMethod`,

and the following functions

* `nextevent_bps`,
* `nextevent_zz`.

The `NextEvent` type encapsulates an object returned when sampling from the IPP is required. It contains:

* a bouncing time
* a function returning whether the bouncing time should be accepted or not (see `Thinning`)
* a flipindex (see `nextevent_zz`)

The functions `nextevent_*` are overloaded for the different possible sampling cases.

### Exact sampling

Exact sampling is possible for specific distributions. In that case, the `dobounce` function returns `true` all the time.
An example is the `MvGaussian` model for which we can indeed sample from the corresponding IPP analytically. See for example the definition of

```julia
nextevent_bps{T<:Vector{Float}}(g::MvGaussian, x::T, v::T)::NextEvent
```

Note the signature of the function. The first parameter indicates how the sampling should be done and the information to do so.
The second and third parameters indicate the ray along which the bouncing time should be produced.

### Sampling via thinning

The `LinearBound` type allows to define a linear upper bound on the intensity of the IPP when you have access to a global bound on the eigenvalues of the Hessian (see [this paper](https://arxiv.org/pdf/1701.04244.pdf) for more details). An accept reject step can then be performed in the `nextevent_*` function (`dobounce`).

An example is the Bayesian logistic regression for which you do have such a bound.

In that case the `nextevent_bps` returns a bouncing time computed thanks to the upper bound and the `dobounce` corresponds to a simple accept/reject step.

## Kernels

The `kernels.jl` code exports a few functions, all of which define how the velocity should be updated in a variety of circumstances. Some of these functions have an exclamation mark attached to them indicating that the transformation is done *in place* (for efficiency). The `refresh_*` functions help indicate how the velocity should be refreshed (e.g. drawn from a spherical gaussian). The `reflect_*` indicate how the velocity should bounce (e.g. specular reflection).
