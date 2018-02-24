export
    Event,
    EventList,
    AllEventList,
    getevent,
    getlastevent,
    pushevent!,
    samplelocalpath,
    quadpathpoly,
    pathmean,
    tmax

"""
    AllowedVarType

Syntactic sugar to make it explicit what variable types are handled in the
local BPS (i.e, the type of the variables on each node of the factor graph).
See also, Event, EventList.
"""
const AllowedVarType = Union{Real, Vector{<:Real}}

"""
    Event{T <: AllowedVarType}

Encapsulation of a local event: a tuple made out of (position, velocity, time)
for a node in the local BPS algorithm.
See also: `EventList`.
"""
struct Event{T <: AllowedVarType}
    x::T
    v::T
    t::Real
end

"""
    EventList{T <: AllowedVarType}

Storage for all events that happened for a variable in the local BPS.
Each variable is stored explicitly (as opposed to storing a list of Event)
so that the times can be traversed efficiently. The xs, vs must all have the
same type within the list.
See also: `Event`, `AllEventList`.
"""
type EventList{T <: AllowedVarType}
    xs::Vector{T}     # coordinates
    vs::Vector{T}     # associated velocities
    ts::Vector{Real}  # associated times
end

function EventList(T::Type)
    @assert T <: AllowedVarType
    EventList(Vector{T}(0), Vector{T}(0), Vector{Real}(0))
end

"""
    AllEventList

Storage for all the `EventList` objects corresponding to a factor graph in a
local BPS setting. The element `evl` can be accessed by the variable index so
`a.evl[k]` will return the eventlist associated with variable `k`. Note that
the different `EventList` may have different types (in particular, different
dimensions). See also: `Event`, `EventList`.
"""
mutable struct AllEventList
    evl::Vector{EventList}
    types::Vector{DataType}
    function AllEventList(types::Vector{DataType})
        new([EventList(types[i]) for i ∈ 1:length(types)], types)
    end
end
AllEventList(T::DataType, nvars::Int) = AllEventList([T for i ∈ 1:nvars])

tmax(aev::AllEventList) = maximum(aev.evl[i].ts[end] for i ∈ 1:length(aev.evl))

"""
    getevent(evl, evidx)

Retrieve a specfic event (at index `evidx`) in an `EventList` object `evl` and
return it as an `Event` object.
See also: `Event`, `EventList`.
"""
function getevent(evl::EventList, evidx::Int)
    Event(evl.xs[evidx], evl.vs[evidx], evl.ts[evidx])
end

"""
    getlastevent(evl)

Retrieve the last `Event` in an `EventList` object `evl`.
See also: `getevent`, `EventList`.
"""
function getlastevent(evl::EventList)
    Event(evl.xs[end], evl.vs[end], evl.ts[end])
end

"""
    pushevent!(evl, ev)

Append an `Event` object `ev` to an `EventList` object `evl`. The `Event`
object is unpacked and pushed to the vectors encapsulated by `evl`.
See also: `EventList`.
"""
function pushevent!(evl::EventList, ev::Event)
    push!(evl.xs, ev.x)
    push!(evl.vs, ev.v)
    push!(evl.ts, ev.t)
    evl
end

"""
    getlocalsegment(evl, j)

Retrieve the `j`th event in an `EventList` object `evl` and return as a tuple
that event and the next event. If there is no next event (end of event list),
an "infinity event" is returned (in order to allow extrapolation).
See also: `samplelocalpath`.
"""
function getlocalsegment(evl::EventList, j::Int)
    if j == length(evl.ts)
        return (getevent(evl,j), Event(NaN, NaN, Inf))
    end
    (getevent(evl,j), getevent(evl,j+1))
end

"""
    samplelocalpath(evl, t)

Sample x_k(t) given the corresponding `EventList` object `evl`. The times `t`
must be a list of sorted times (or a single positive time).
See also: `getlocalsegment`.
"""
function samplelocalpath(evl::EventList, t::AllowedTimeType)
    @assert issorted(t) "list of times should be sorted"
    ti  = t[1]
    @assert ti >= 0.0 "the first time should be greater than 0.0"
    # find the starting local segment
    j = searchsortedlast(evl.ts, ti)
    # event and next event (describes the "segment")
    eva, evb = getlocalsegment(evl,j)
    # allocation of samples
    samples = Vector{typeof(eva.x)}(length(t))
    # recursively work with local segments
    i=1
    while i<= length(t)
        ti = t[i]
        if ti <= evb.t
            # stay in the current segment
            samples[i] = eva.x + (ti - eva.t) * eva.v
            i += 1
        else
            # get the next segment
            j += searchsortedlast(evl.ts[j+1:end], ti)
            eva, evb = getlocalsegment(evl, j)
        end
    end
    samples
end
samplelocalpath(evl::EventList, t::Real) = samplelocalpath(evl, [t])[1]

"""
    quadpathpoly(evl, pol, T)

Integrate a polynomial on the elements along a path for a single node.
So if we are sampling in Rd and have phi(X)=Poly(xk), then we can integrate the
polynomial along the path. Note that we can't do the same thing if phi(X) mixes
2 or more nodes. See also syntactic sugar to perform the same operation on all
nodes.
"""
function quadpathpoly(evl::EventList, pol::Poly, T::Real)
    dim  = length(evl.xs[1])
    res  = zeros(dim)
    nseg = searchsortedlast(evl.ts,T)
    for segidx = 1:nseg-1
        # for each segment in the path
        eva, evb = getlocalsegment(evl, segidx)
        # segment time, initial point and velocity
        tau = evb.t - eva.t
        xa  = eva.x
        v   = (evb.x - xa) / tau
        # integrate along that segment
        for d = 1:dim
            res[d] += polyint(polyval(pol,Poly([xa[d],v[d]])))(tau)
        end
    end
    # last segment
    eva = getevent(evl, nseg)
    # characteristics
    tau = T - eva.t
    for d = 1:dim
        res[d] += polyint(polyval(pol,Poly([eva.x[d],eva.v[d]])))(tau)
    end
    return res / T
end
function quadpathpoly(aev::AllEventList, pol::Poly, T::Real)
    res = Vector{Vector{<:Real}}(length(aev.evl))
    for (d, evl) in enumerate(aev.evl)
        res[d] = quadpathpoly(evl, pol, T)
    end
    return res
 end
pathmean(evl::EventList, T::Real) = quadpathpoly(evl, Poly([0.0,1.0]), T)
pathmean(aev::AllEventList, T::Real) = quadpathpoly(aev, Poly([0.0,1.0]), T)
pathmean(aev::AllEventList) = pathmean(aev, tmax(aev))
