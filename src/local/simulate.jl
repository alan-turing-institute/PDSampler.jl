export
    LocalSimulation,
    simulate
"""
    LocalSimulation

Describe a Local-BPS sampling
"""
immutable LocalSimulation
    fg::FactorGraph
    x0::Vector{AllowedVarType}  # Starting point
    v0::Vector{AllowedVarType}  # Starting velocity
    T::AbstractFloat                    # Simulation time
    maxnevents::Int             # Maximum number of events
    lambdaref::AbstractFloat            # Refreshment rate
    # -- implicit
    nvars::Int                  # Number of nodes
    nfactors::Int               # Number of factors
    function LocalSimulation(fg, x0, v0, T, maxnevents, lambdaref)
        # check none of the default named arguments went through
        cond =  !(  fg == :undefined ||
                    x0 == :undefined ||
                    v0 == :undefined ||
                    T  == :undefined )
        @assert cond "One or more essential argument undefined"
        cond = (length(x0) == length(v0) == fg.structure.nvars > 0)
        @assert cond "inconsistent dims"
        cond = (T > 0.0 && maxnevents > 0 && lambdaref > 0.0)
        @assert cond "inconsistent params"
        # creating the object
        new( fg, x0, v0, T,
             maxnevents,
             lambdaref,
             length(x0),
             length(fg.factors) )
    end
end
# Constructor with named arguments
function LocalSimulation(;
            factorgraph = :undefined,
            x0          = :undefined,
            v0          = :undefined,
            T           = :undefined,
            maxnevents  = 50000,
            lambdaref   = 0.5)::LocalSimulation
    # calling the unnamed constructor
    LocalSimulation( factorgraph,
                     x0, v0, T,
                     maxnevents,
                     lambdaref )
end

"""
    simulate(sim)

Run a local BPS simulation following the specifications of `sim`.
"""
function simulate(sim::LocalSimulation)::Tuple{AllEventList, Dict}

    (start, all_evlist, pq, tref) = ls_init(sim)

    # counters for the first branch (bounce) and the second branch (refresh)
    ev_firstbranch  = 0
    ev_secondbranch = 0
    nevents         = 0 # sum of above two
    # global clock
    globalclock     = 0.0

    # global iteration
    prog = Progress(sim.maxnevents, 1)
    # HACK
    #while globalclock < sim.T && nevents < sim.maxnevents
    for evnum in 1:sim.maxnevents
        # get bounce and dequeue
        (fidx, tbounce) = peek(pq)
        dequeue!(pq)
        t = tbounce
        globalclock += t
        # ------------
        # FIRST BRANCH (t<tref => )
        if t < tref
            ev_firstbranch += 1
            ls_firstbranch!(sim.fg, fidx, all_evlist, pq, t)
        # -------------
        # SECOND BRANCH (t>=tref => refreshment)
        else
            ev_secondbranch += 1
            pq    = ls_refreshment(sim.fg, tref, all_evlist)
            tref += randexp()/sim.lambdaref
        end
        # if use while loop
        # nevents += 1
        next!(prog)
    end

    details = Dict(
        "clocktime"       => time()-start,
        "nevents"         => nevents,
        "ev_firstbranch"  => ev_firstbranch,
        "ev_secondbranch" => ev_secondbranch,
        "globalclock"     => globalclock
    )

    (all_evlist, details)
end

# =============================================================================
# NOTE: The functions with a `ls_` prefix are helper functions exclusively
# designed for the local simulation. They aren't expected to be used elsewhere.

"""
    ls_init(sim)

Initialise a local simulation: store the starting clock time, create a list of
all eventlists and push the original event, create the priorityqueue and fill
it with the initial bounce times
"""
function ls_init(sim::LocalSimulation
                )::Tuple{AbstractFloat,AllEventList,PriorityQueue,AbstractFloat}
    # instantiate wall clock
    start = time()
    # initialisation of the eventlists for every node in the graph
    all_evlist = AllEventList([typeof(sim.x0[i]) for i in 1:sim.nvars])
    for i in 1:sim.nvars
        evi = Event(sim.x0[i], sim.v0[i], 0.0)
        pushevent!(all_evlist.evl[i], evi)
    end
    # initialisation of the priority queue and the refreshment time
    pq   = PriorityQueue(Int, AbstractFloat)
    tref = randexp()/sim.lambdaref
    # filling of the priority queue with initial position
    for fidx in 1:sim.fg.structure.nfactors
        (xf, vf, g) = ls_retrieve(sim.fg, fidx, all_evlist, 0.0)
        ls_updatepq!(pq, sim.fg, fidx, xf, vf, g, 0.0)
    end
    (start, all_evlist, pq, tref)
end

"""
    ls_firstbranch(fg, fidx, all_evlist, pq, t)

Operation corresponding to the first branch of events in simulate.
"""
function ls_firstbranch!(fg::FactorGraph, fidx::Int, all_evlist::AllEventList,
                         pq::PriorityQueue, t::AbstractFloat
                         )::Tuple{AllEventList,PriorityQueue}
    # retrieve xf, vf corresponding to factor
    (xf, vf, g, vars) = ls_retrieve(fg, fidx, all_evlist, t, true)
    ls_saveupdate!(all_evlist, vars, xf, vf, t)
    ls_updatepq!(pq, fg, fidx, xf, vf, g, t)
    # same story for linked factors (fp)
    for fpidx in linkedfactors(fg, fidx)
        # we don't need to retrieve `vars` here
        (xfp, vfp, gp) = ls_retrieve(fg, fpidx, all_evlist, t)
        ls_updatepq!(pq, fg, fpidx, xfp, vfp, gp, t)
    end
    (all_evlist, pq)
end


"""
    ls_reshape(u, v)

(Local Simulation, helper function) Reshape a single vector `u` to a vector of
AllowedVarType following the structure of `v`. For example, let
`v = [0., [0., 0.]]` and `u=[1., 2., 3.]` then `ls_reshape(u,v)` will return
`[1., [2., 3.]]`.
"""
function ls_reshape{V<:Vector{AllowedVarType}}(u::Vector{AbstractFloat}, v::V)::V
    w = similar(v)
    # offsets to know where to look for next block of information
    offset = 0
    for i in 1:length(v)
        tmpi = length(v[i])
        # if length is 1, return only single point, otherwise vector
        @inbounds w[i] = (tmpi==1) ? u[offset+1] : u[offset+(1:tmpi)]
        offset += tmpi
    end
    w
end

"""
    ls_retrieve(fg, fidx, all_evlist, t, doreflect)

(Local simulation, helper function) Retrieve the structure corresponding to all
positions of the nodes attached to a factor (`xf`) and the velocities (`vf`) as
well as the gradient at `xf` and the list of indices of the attached variables.
"""
function ls_retrieve(fg::FactorGraph, fidx::Int,
                    all_evlist::AllEventList,
                    t::AbstractFloat, doreflect::Bool=false)
    # indices of the variable associated with factor fidx
    vars = assocvariables(fg, fidx)
    # allocate xf, vf, note the different variables
    # don't necessarily have the same types
    vf = Vector{AllowedVarType}(length(vars))
    xf = similar(vf)
    # shortcut to initial event
    if t == 0.0
        for (i, k) in enumerate(vars)
            ev = getevent(all_evlist.evl[k], 1)
            xf[i], vf[i] = ev.x, ev.v
        end
    # retrieving later than initial event (standard case)
    else
        # retrieve vf, xf(t)
        for (i, k) in enumerate(vars)
            # get the eventlist corresponding to variable k
            ev    = getlastevent(all_evlist.evl[k])
            # store the components and extrapolate from ev.x as needed
            vf[i] = ev.v
            xf[i] = ev.x + (t-ev.t)*ev.v
        end
    end
    # compute the gradient of the factor's log likelihood at xf
    g  = fg.factors[fidx].gll(vcat(xf...))
    # compute the reflection (if required)
    vf = doreflect ? ls_reshape(reflect_bps!(g, vcat(vf...)), vf) : vf
    (xf, vf, g, vars)
end

"""
    ls_saveupdate!(all_evlist, vars, xf, vf, t)

(Local simulation, helper function) Append the information contained in `xf`,
`vf` in the appropriate EventLists corresponding to the relevant variables
indexed in `vars` (variables associated with the factor `f`, `xf` and `vf` are
the recovered positions and velocities obtained when considering factor `f`).
"""
function ls_saveupdate!(all_evlist::AllEventList, vars::Vector{Int},
                        xf::Vector{AllowedVarType}, vf::Vector{AllowedVarType},
                        t::AbstractFloat)::AllEventList
    # Add the new sample to the list
    for (i, k) in enumerate(vars)
        pushevent!(all_evlist.evl[k], Event(xf[i], vf[i], t))
    end
    all_evlist
end

"""
    ls_updatepq!(pq, fg, fidx, xf, vf, g, t)

(Local simulation, helper function) Compute the bounce time (possibly via
thinning depending on the factor's implementation) associated with a factor and
push it to the PriorityQueue `pq`.
"""
function ls_updatepq!(pq::PriorityQueue, fg::FactorGraph, fidx::Int,
                      xf::Vector{AllowedVarType}, vf::Vector{AllowedVarType},
                      g::Vector{AbstractFloat}, t::AbstractFloat)::PriorityQueue
    # useful temporary variables
    vcxf, vcvf = vcat(xf...), vcat(vf...)
    # Update time in Priority Queue for the current factor
    acc  = false
    tauf = 0.0
    while !acc
        bounce    = fg.factors[fidx].nextevent(vcxf, vcvf)
        acc, tauf = bounce.dobounce(g, vcvf), bounce.tau
    end
    # add the time to current priorityQueue
    pq[fidx] = t+tauf
    pq
end

"""
    ls_random(evl)

(Local simulation, helper function) Generate N(0,1) random numbers of dimension
corresponding to the passed `EventList`.
"""
function ls_random(evl::EventList)::AllowedVarType
    l = length(evl.xs[end])
    l>1?randn(l):randn()
end

"""
    ls_refreshment(fg, t, all_evlist)

(Local simulation, helper function) Generate a new PriorityQueue by drawing
new velocities for each node and storing it (cf. NewQueue algorithm in the BPS
paper).
"""
function ls_refreshment(fg::FactorGraph, t::AbstractFloat,
                        all_evlist::AllEventList)::PriorityQueue
    # draw a new bunch of velocities for each node
    v = Vector{AllowedVarType}(fg.structure.nvars)
    for i in 1:fg.structure.nvars
        @inbounds v[i] = ls_random(all_evlist.evl[i])
    end
    # Instantiate a new priority queue
    pq = PriorityQueue(Int, AbstractFloat)
    for fidx in 1:fg.structure.nfactors
        # retrieve xf, vf corresponding to factor
        (xf, vf, g, vars) = ls_retrieve(fg, fidx, all_evlist, t)
        # update the priority queue, note v[vars] and not vf
        # since we are using the refreshed velocity
        ls_updatepq!(pq, fg, fidx, xf, v[vars], g, t)
    end
    pq
end
