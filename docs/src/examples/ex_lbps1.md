# Local BPS (Chain of Gaussians)

(*the code for this example can be found [here](https://github.com/alan-turing-institute/PDSampler.jl/blob/master/test/ex_lbps1.jl), note that the doc rendered here was automatically generated, if you want to fix it, please do it in the julia code directly*)

The approach to using the local BPS is much the same as for the global one
except that you need to specify a `FactorGraph`.
That object will contain the structure of the factor graph (i.e.: which factor is connected to which variables) as well as the list of all factors (which have a `lgradll` and `nextevent` since each factor can be seen individually as a small BPS).

Below, we show how to declare a chain of bivariate gaussians:

```julia
using PDSampler, Random
nfac = 3 # number of factors

mvg = MvGaussianStandard(zeros(2),diagm(0=>ones(2)))

# all factors have that same likelihood
chainfactor(i) = Factor(
                    (x,v)->nextevent_bps(mvg, x, v),
                    (x)->gradloglik(mvg, x),
                    i )

# assemble into a chain graph
chain = chaingraph([chainfactor(i) for i in 1:nfac])

```
This is a simple graph with a known structure so that it's already defined
through the `chaingraph` function (in `src/local/factorgraph.jl`). For an
arbitrary graph, you would need to provide two things:
1. the structure of the factor graph: a list of list where each element corresponds to a factor and the corresponding list contains the indices of the variables attached to that factor
2. the list of factors
The rest is very similar to the global BPS:

```julia
Random.seed!(123)
lambdaref  = .01
maxnevents = 10000
T          = Inf
nvars      = chain.structure.nvars
x0         = randn(nvars)
v0         = randn(nvars)
v0        /= norm(v0)

lsim = LocalSimulation(chain, x0, v0, T, maxnevents, lambdaref)

(all_evlist, details) = simulate(lsim)

```
The `all_evlist` object contains a list of `EventList` corresponding to what
happened on each of the factors. It can also be sampled using `samplelocalpath`
(cf. `src/local/event.jl`).
