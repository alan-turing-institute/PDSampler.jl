# [Core tools](@id td-coretools)

## Geometry

The `geometry.jl` code exports the following **immutable** types:

* `Unconstrained <: Domain `,
* `Polygonal <: Domain`.

and a function

* `nextboundary`

Any geometry type needs to have a `nextboundary` function associated with it of the form `nextboundary(g, x, v)` where `g` refers to the geometry, `x` to the position and `v` to the velocity. The aim of that function is to:

* return a time of first hit $\tau_h$ (following the ray $x+vt$, for $t>0$),
* return the normal to the boundary at that hitting point.

So for example, the `Unconstrained` is an empty type and the associated `nextboundary` function just returns `(NaN, NaN)` indeed a boundary will never be crossed and there is no normal. These `NaN` are processed by calling functions.

The `Polygonal` domain requires the definition of the normals and the intercepts. The `nextboundary` function is pretty simple (intersection of lines). Note a few tricks for numerical stability:

* near parallel case can be ignored (the crossing time will be huge compared to bouncing time or refreshment time and therefore we can just ignore it)
* negative times can be ignored (we're only going forward)
* times that are very close to zero are ignored (means that we are currently already very close to the boundary meaning that we will bounce away)

## Sampling from an IPP

The `ippsampler.jl` exports the following **immutable** types:

* `NextEvent`
* `LinearBound <: Thinning <: IPPSamplingMethod` 

## Kernels
