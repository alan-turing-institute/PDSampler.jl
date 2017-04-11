# PDMP Types

(**/!\WIP/!\**)

```@docs
FactorGraph
```

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
