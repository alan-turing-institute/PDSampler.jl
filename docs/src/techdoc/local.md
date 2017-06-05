# [Local sampler](@id td-localsampler)

Link to the source files:

* [`local/event.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/local/event.jl)
* [`local/factorgraph.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/local/factorgraph.jl)
* [`local/simulate.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/local/simulate.jl)

## Event

### Hierarchy of types

The events in the local settings are triples of the form `(x,v,t)`. These can be encapsulated in the immutable `Event` type.
Every node in the graph has a list of such events attached to it.
For efficiency reasons, the list of event is not actually a list of `Event` but rather another type `EventList` which allows traversing the corresponding data structure efficiently.

The list of `EventList` (i.e. all events that have happened on the graph) is encapsulated in another type `AllEventList`.


#### EventList

The `EventList` type corresponds to what is stored by a single node of the factor graph. It contains the list of positions `xs`, the list of velocities `vs` and the list of times `ts` associated with that node.

The function `getevent`, `getlastevent`, `pushevent!`, `getlocalsegment`, `samplelocalpath`, `quadpathpoly` and `pathmean` all work on an `EventList` object.

In other words, an `EventList` is analogous to the `Path` object of the global sampler.

#### AllEventList

This is just a vector of `EventList`. It also keeps track of the types associated with each `EventList`. Indeed, each node may correspond to a variable of different dimensionality so these `types` make initialisation procedures simpler.

### Auxiliary functions

* `getevent` retrieve an event in an `EventList`
* `getlastevent` retrieve the last event of an `EventList`
* `pushevent!` add an event to an `EventList`
* `getlocalsegment` get an event and the subsequent event
* `samplelocalpath` sample the variable corresponding to a node by taking samples along the path described by the corresponding `EventList`
* `quadpathpoly` integrate a polynomial along the path described by an `EventList`

## Factor Graph

### Hierarchy of types

All types are immutables (define the graph).

* The base type is a `Factor`
* The connection pattern is a `FactorGraphStruct`
* The list of factors + structure forms a `FactorGraph`

#### Factor

A `Factor` encapsulates

* `nextevent` a function which is able to produce a first arrival time from the IPP corresponding to that factor
* `gll` the gradient of the loglikelihood attached to that factor
* `index` this is a dummy index to be able to refer to specific factors in the factor graph

#### FactorGraphStruct

A `FactorGraphStruct` encapsulates

* `flist` a list of list, every entry corresponds to a factor and the list of variables that factor is connected to.

implicitly it also encapsulates

* `vlist` a list of list, every entry corresponds to a variable and the list of factors that variable is connected to.
* `nfactors`, `nvars` the number of factors and variables

The function `chainstruct` defines a returns a `FactorGraphStruct` corresponding to a chain with a given number of nodes (variables).

#### FactorGraph

A `FactorGraph` encapsulates

* `structure` a `FactorGraphStruct` object
* `factors` a vector of `Factor` objects

### Auxiliary functions

* `assocvariables` and `assocfactors` return the indices of the associated variables to a factor (resp. factors to a variable)
* `linkedfactors` for a given factor returns the list of factors which share a variable with it

## Simulate

That file has the same structure as that for the global sampler with one major `LocalSimulation` object encapsulating all the parameters of a given simulation.

The main loop uses a number of helper functions `ls_*` in order to make the logic appear more clearly.
Focusing on the main loop for now, it should be clear that is centered around a priority queue.

```julia
(fidx, tbounce) = peek(pq)
...
if t < tref
    ...
else
    ...
end
```

So the minimum time is recovered from the priority queue as well as the index of the corresponding factor. Then, there is a check to see whether we are in the bouncing scenario (*first branch*) or the refreshment scenario (*second branch*).

### Helper functions

#### Initialisation

The function `ls_init` initialises an `AllEventList` object corresponding to the graph as well as a priority queue.

#### Reshape

Once a factor is selected, the operations are very similar to the global BPS. But once an event is generated, it is necessary to slice it according to the different nodes so that they all store the relevant portion of information. The `reshape` function reshapes a generated object into the appropriate format.

#### Retrieve

The retrieve function is called to access a specific factor. For that factor it retrieves:

* The positions `xf`, a list of positions for each nodes associated with the factor interpolated at a time `t`
* The velocities `vf`, a list of velocities for each nodes associated with the factor
* The gradient at `xf` (as seen from the factor)
* The index of the variables associated with the factor

The interpolation is done following the method highlighted in the paper (recover the last event and follow the ray for required time). 
