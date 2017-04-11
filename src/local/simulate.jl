using Base.Collections:
        PriorityQueue,
        enqueue!,
        dequeue!,
        peek

export
    LocalSimulation,
    simulate,
    # helper functions
    ls_reshape,
    ls_retrieve,
    ls_saveupdate!,
    ls_updatepq!,
    ls_random,
    ls_refreshment

"""
    LocalSimulation

Describe a Local-BPS sampling
"""
immutable LocalSimulation
    fg::FactorGraph
    x0::Vector{AllowedVarType}  # Starting point
    v0::Vector{AllowedVarType}  # Starting velocity
    T::Float                    # Simulation time
    maxnevents::Int             # Maximum number of events
    lambdaref::Float            # Refreshment rate
    # -- implicit
    nvars::Int                  # Number of nodes
    function LocalSimulation(fg, x0, v0, T, maxnevents, lambdaref)
        @assert length(v0)==length(x0)==fg.structure.nvars "inconsistent dims"
        @assert T>0.0 && maxnevents > 0 && lambdaref > 0.0 "inconsistent params"
        new(fg,x0,v0,T,maxnevents,lambdaref, length(x0))
    end
end

function simulate(sim::LocalSimulation)::Tuple{AllEventList, Dict}
    # keep track of how much time we've been going for
    start = time()

    # initialisation of the eventlists for every node in the graph
    all_evlist = AllEventList([typeof(sim.x0[i]) for i in 1:sim.nvars])

    for i in 1:sim.nvars
        evi = Event(sim.x0[i], sim.v0[i], 0.0)
        pushevent!(all_evlist.evl[i], evi)
    end

    # initialisation of the priority queue and the refreshment time
    pq   = PriorityQueue(Int, Float)
    tref = randexp()/sim.lambdaref
    # filling of the priority queue with initial position
    for (fidx, factor) in enumerate(sim.fg.factors)
        vars   = assocvariables(sim.fg, fidx)
        xf, vf = sim.x0[vars], sim.v0[vars]
        g      = factor.gll(vcat(xf...))
        pq     = ls_updatepq!(pq, sim.fg, fidx, xf, vf, g, 0.0)
    end

    # counters for the first branch (bounce) and the second branch (refresh)
    ev_firstbranch  = 0
    ev_secondbranch = 0
    nevents         = 0 # sum of above two
    # global clock
    globalclock = 0.0

    # global iteration
    while globalclock < sim.T && nevents < sim.maxnevents
        # get bounce and dequeue
        (fidx, tbounce) = peek(pq)
        dequeue!(pq)
        t = tbounce
        globalclock += t
        # ------------
        # FIRST BRANCH (t<tref => )
        if t < tref
            ev_firstbranch += 1
            # retrieve xf, vf corresponding to factor
            (xf, vf, g, vars) = ls_retrieve(sim.fg, fidx, all_evlist, t, true)
            ls_saveupdate!(all_evlist, vars, xf, vf, t)
            ls_updatepq!(pq, sim.fg, fidx, xf, vf, g, t)
            # same story for linked factors (fp)
            for fpidx in linkedfactors(sim.fg, fidx)
                (xfp, vfp, gp) = ls_retrieve(sim.fg, fpidx, all_evlist, t)
                ls_updatepq!(pq, sim.fg, fpidx, xfp, vfp, gp, t)
            end
        # -------------
        # SECOND BRANCH (t>=tref => refreshment)
        else
            ev_secondbranch += 1
            pq    = ls_refreshment(sim.fg, tref, all_evlist)
            tref += randexp()/sim.lambdaref
        end
        nevents += 1
    # end of nevents = 1:sim.maxnevents
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

# ==============================================================================
# NOTE: The functions with a `ls_` prefix are helper functions exclusively
# designed for the local simulation.

"""
    ls_reshape(u, v)

(Local Simulation, helper function) Reshape a single vector `u` to a vector of
AllowedVarType following the structure of `v`. For example, let
`v = [0., [0., 0.]]` and `u=[1., 2., 3.]` then `ls_reshape(u,v)` will return
`[1., [2., 3.]]`.
"""
function ls_reshape{V<:Vector{AllowedVarType}}(u::Vector{Float}, v::V)::V
    w = Vector{AllowedVarType}(length(v))
    offset = 0
    for i in 1:length(v)
        tmpi    = length(v[i])
        # if length is 1, return only single point, otherwise vector
        w[i]    = (tmpi==1) ? u[offset+1] : u[offset+(1:tmpi)]
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
                    t::Float, doreflect::Bool=false)
    # indices of the variable associated with factor fidx
    vars  = assocvariables(fg, fidx)
    # allocate xf, vf, note the different variables
    # don't necessarily have the same types
    vf = Vector{AllowedVarType}(length(vars))
    xf = Vector{AllowedVarType}(length(vars))
    # retrieve vf, xf(t)
    for (i, k) in enumerate(vars)
        # get the eventlist corresponding to variable k
        ev    = getlastevent(all_evlist.evl[k])
        # store the components and extrapolate from ev.x as needed
        vf[i] = ev.v
        xf[i] = ev.x + (t-ev.t)*ev.v
    end
    # compute the gradient of the factor's log likelihood at xf
    g  = fg.factors[fidx].gll(vcat(xf...))
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
                        t::Float)::AllEventList
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
                      g::Vector{Float}, t::Float)::PriorityQueue
    # useful temporary variables
    vcxf, vcvf = vcat(xf...), vcat(vf...)
    # Update time in Priority Queue for the current factor
    bounce = fg.factors[fidx].nextevent(vcxf, vcvf)
    acc    = bounce.dobounce(g, vcvf)
    while !acc
        bounce = fg.factors[fidx].nextevent(vcxf, vcvf)
        acc    = bounce.dobounce(g, vcvf)
    end
    tauf = bounce.tau
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
function ls_refreshment(fg::FactorGraph, t::Float,
                        all_evlist::AllEventList)::PriorityQueue
    # draw a new bunch of velocities for each node
    v = Vector{AllowedVarType}(fg.structure.nvars)
    for i in 1:fg.structure.nvars
        v[i] = ls_random(all_evlist.evl[i])
    end
    # Instantiate a new priority queue
    pq = PriorityQueue(Int,Float)
    for fidx in 1:fg.structure.nfactors
        # retrieve xf, vf corresponding to factor
        (xf, vf, g, vars) = ls_retrieve(fg, fidx, all_evlist, t)
        # update the priority queue, note v[vars] and not vf
        # since we are using the refreshed velocity
        ls_updatepq!(pq, fg, fidx, xf, v[vars], g, t)
    end
    pq
end
