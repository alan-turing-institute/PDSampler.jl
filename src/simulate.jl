export
    Simulation,
    simulate

"""
    Simulation

Describe a PDMP sampling: information about the initial point, the time of the
simulation, the function to sample from an IPP, etc.
"""
immutable Simulation
    x0::Vector{Float}      # Starting point
    v0::Vector{Float}      # Starting velocity
    T::Float               # Simulation time
    nextevent::Function    # Appropriate simulation for first arrival time
    gll::Function          # Gradient of Log Lik (potentially CV)
    nextboundary::Function # Where/When is the next boundary hit
    lambdaref::Float       # Refreshment rate
    algname::String        # BPS, ZZ

    # derived
    dim::Int             # dimensionality

    # optional named arguments
    mass::Matrix{Float}  # mass matrix (preconditioner)
    blocksize::Int       # increment the storage by blocks
    maxsimtime::Float    # max. simulation time (s)
    maxsegments::Int     # max. num. of segments
    maxgradeval::Int     # max. num. grad. evals
    refresh!::Function   # refreshment function (TODO: examples)

    # constructor
    function Simulation( x0, v0, T, nextevent, gradloglik, nextboundary,
                lambdaref = 1.0, algname = "BPS";
                mass = eye(0), blocksize = 1000, maxsimtime = 4e3,
                maxsegments = Int(1e6), maxgradeval = Int(1e8),
                refresh! = refresh_global! )
        #
        an = uppercase(algname)
        @assert (an in ["BPS", "ZZ"]) "Unknown algorithm <$algname>"
        new( x0, v0, T, nextevent, gradloglik, nextboundary, lambdaref,
             an, length(x0), mass, blocksize, maxsimtime, maxsegments,
             maxgradeval, refresh! )
    end
end

"""
    simulate(sim)

Launch a PDMP simulation defined in the sim variable. Return the corresponding
Path and a dictionary of indicators (clocktime, ...).
"""
function simulate(sim::Simulation)::Tuple{Path, Dict}
    # keep track of how much time we've been going for
    start = time()

    # dimensionality
    d = sim.dim
    # initial states
    x, v = copy(sim.x0), copy(sim.v0)
    # time counter, and segment counter
    t, i = 0.0, 1
    # counters for the number of effective loops and
    # for the number of evaluations of the gradient
    # this will be higher than the number of segments i
    lcnt, gradeval = 0, 0

    # storing by blocks of blocksize nodes at the time
    blocksize = sim.blocksize

    # storing xs as a single column for efficient resizing
    xs, ts = zeros(d*blocksize), zeros(blocksize)
    # store initial point
    xs[1:d] = x
    # mass matrix?
    mass = copy(sim.mass) # store it here as we may want to adapt it

    # compute current reference bounce time
    lambdaref = sim.lambdaref
    tauref    = randexp()/lambdaref
    # Compute time to next boundary + normal
    (taubd, normalbd) = sim.nextboundary(x, v)

    # keep track of how many refresh events
    nrefresh = 0
    # keep track of how many boundary hits
    nboundary = 0
    # keep track of how many standard bounces
    nbounce = 0

    while t < sim.T && gradeval < sim.maxgradeval
        # increment the counter to keep track of the number of effective loops
        lcnt += 1
        # simulate first arrival from IPP
        bounce = sim.nextevent(x, v)
        # find next event (unconstrained case taubd=NaN (ignored))
        tau = min(bounce.tau,taubd,tauref)

        # standard bounce
        if tau == bounce.tau
            # there will be an evaluation of the gradient
            gradeval += 1
            nbounce  += 1
            # updating position/time
            t += tau
            x += tau*v
            # exploiting the memoryless property
            tauref -= tau
            # ---- BOUNCE ----
            g = sim.gll(x)
            if bounce.dobounce(g, v) # e.g.: thinning, acc/rej
                if sim.algname == "BPS"
                    if length(mass)>0
                        v = reflect_bps!(g, v, mass)
                    else
                        v = reflect_bps!(g,v)
                    end
                elseif sim.algname == "ZZ"
                    v = reflect_zz!(bounce.flipindex, v)
                end
            else
                # move closer to the boundary/refreshment time
                taubd  -= tau
                # we don't need to record when rejecting
                continue
            end
        # hard bounce against boundary
        elseif tau == taubd
            #
            nboundary += 1
            # Record point epsilon from boundary for numerical stability
            x += (tau-1e-10)*v
            t += tau
            # exploiting the memoryless property
            tauref -= tau
            # ---- BOUNCE (boundary) ----
            if sim.algname == "BPS"
                if length(mass)>0
                    v = reflect_bps!(normalbd, v, mass)
                else
                    v = reflect_bps!(normalbd,v)
                end
            elseif sim.algname == "ZZ"
                v = reflect_zz!(find((v.*normalbd).<0.0), v)
            end
        # random refreshment
        else
            #
            nrefresh += 1
            #
            x += tau*v
            t += tau
            # ---- REFRESH ----
            v = sim.refresh!(v)
            # update tauref
            tauref = randexp()/lambdaref
        end
        # check when the next boundary hit will occur
        (taubd, normalbd) = sim.nextboundary(x, v)
        # increment the counter for the number of segments
        i += 1

        # Increase storage on a per-need basis.
        if mod(i,blocksize)==0
            resize!(xs,length(xs)+d*blocksize)
            resize!(ts,length(ts)+blocksize)
        end

        # Storing path times
        ts[i] = t
         # storing columns vertically, a vector is simpler/cheaper to resize
        xs[((i-1)*d+1):(i*d)] = x

        # Safety checks to break long loops every 100 iterations
        if mod(lcnt,100)==0
            if  (time()-start) > sim.maxsimtime
                println("Max simulation time reached. Stopping")
                break
            end
            if  i > sim.maxsegments
                println("Too many segments generated. Stopping")
                break
            end
        end
    end # End of while loop

    details = Dict(
        "clocktime" => time()-start,
        "ngradeval" => gradeval,
        "nloops"    => lcnt,
        "nsegments" => i,
        "nbounce"   => nbounce,
        "nboundary" => nboundary,
        "nrefresh"  => nrefresh
    )

    (Path(reshape(xs[1:(i*d)],(d,i)), ts[1:i]), details)
end
