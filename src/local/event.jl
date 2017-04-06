"""
    AllowedVarType

Syntactic sugar to make it explicit what variable types are handled in the
local BPS (i.e, the type of the variables on each node of the factor graph).
See also, Event, EventList.
"""
AllowedVarType = Union{Float, Vector{Float}}

"""
    Event{T <: AllowedVarType}

Encapsulation of a local event: a tuple made out of (position, velocity, time)
for a node in the local BPS algorithm.
See also: `EventList`.
"""
immutable Event{T <: AllowedVarType}
    x::T
    v::T
    t::Float
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
    ts::Vector{Float} # associated times
    # Constructor
    EventList() = new(Vector{T}(0), Vector{T}(0), Vector{Float}(0))
end

"""
    AllEventList

Storage for all the `EventList` objects corresponding to a factor graph in a
local BPS setting. The element `evl` can be accessed by the variable index so
`a.evl[k]` will return the eventlist associated with variable `k`. Note that
the different `EventList` may have different types (in particular, different dimensions).
See also: `Event`, `EventList`.
"""
type AllEventList
    evl::Vector{EventList}
    types::Vector{DataType}
    function AllEventList(types::Vector{DataType})
        new([EventList{types[i]}() for i in 1:length(types)],
            types)
    end
end
AllEventList(T::DataType, nvars::Int) = AllEventList([T for i in 1:nvars])

"""
    getevent(evl, evidx)

Retrieve a specfic event (at index `evidx`) in an `EventList` object `evl` and
return it as an `Event` object.
See also: `Event`, `EventList`.
"""
function getevent(evl::EventList, evidx::Int)::Event
    Event(evl.xs[evidx],evl.vs[evidx],evl.ts[evidx])
end

"""
    getlastevent(evl)

Retrieve the last `Event` in an `EventList` object `evl`.
See also: `getevent`, `EventList`.
"""
function getlastevent(evl::EventList)::Event
    Event(evl.xs[end],evl.vs[end],evl.ts[end])
end

"""
    pushevent!(evl, ev)

Append an `Event` object `ev` to an `EventList` object `evl`.
See also: `EventList`.
"""
function pushevent!(evl::EventList, ev::Event)::EventList
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
function getlocalsegment(evl::EventList, j::Int)::Tuple{Event,Event}
    if j==length(evl.ts)
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
            samples[i] = eva.x + (ti-eva.t)*eva.v
            i += 1
        else
            # get the next segment
            j += searchsortedlast(evl.ts[j+1:end],ti)
            eva,evb = getlocalsegment(evl,j)
        end
    end
    samples
end
samplelocalpath(evl::EventList, t::Float) = samplelocalpath(evl, [t])[1]
