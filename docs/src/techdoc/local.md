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

## Simulate
