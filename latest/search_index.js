var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#PDMP.jl-Documentation-1",
    "page": "Introduction",
    "title": "PDMP.jl Documentation",
    "category": "section",
    "text": "PDMP.jl is a package designed to provide an efficient, flexible, and expandable framework for samplers based on Piecewise Deterministic Markov Processes and their applications. This includes the Bouncy Particle Sampler and the Zig-Zag Sampler. See the references at the bottom of this page.The package is currently under construction (Spring 2017). The project is hosted and maintained by the Alan Turing Institute (ATI) and currently developed by Thibaut Lienart. If you encounter problems, please open an issue on Github. If you have comments or wish to collaborate, please send an email to tlienart σ turing > ac > uk."
},

{
    "location": "index.html#Using-the-Package-1",
    "page": "Introduction",
    "title": "Using the Package",
    "category": "section",
    "text": "To install the (currently unregistered) package, use the following command inside the Julia REPL:Pkg.clone(\"git://github.com/alan-turing-institute/PDMP.jl.git\")To load the package, use the command:using PDMPYou can also run the tests with Pkg.test(\"PDMP\") and update to the latest Github version with Pkg.update(\"PDMP\")."
},

{
    "location": "index.html#Examples-1",
    "page": "Introduction",
    "title": "Examples",
    "category": "section",
    "text": "The following examples will introduce you to the functionalities of the package.Pages = [\n    \"examples/bps_mvg_constr.md\",\n    \"examples/lbps_gchain.md\"\n    ]\nDepth = 2"
},

{
    "location": "index.html#Code-documentation-1",
    "page": "Introduction",
    "title": "Code documentation",
    "category": "section",
    "text": "These pages introduce you to the core of the package and its interface. This is useful if you are looking into expanding the code yourself to add a capacity or a specific model.Pages = [\n    \"techdoc/types.md\"\n    ]\nDepth = 2"
},

{
    "location": "index.html#References-1",
    "page": "Introduction",
    "title": "References",
    "category": "section",
    "text": "Alexandre Bouchard-Côté, Sebastian J. Vollmer and Arnaud Doucet, The Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain Monte Carlo Method, arXiv preprint, 2015.\nJoris Bierkens, Alexandre Bouchard-Côté, Arnaud Doucet, Andrew B. Duncan, Paul Fearnhead, Gareth Roberts and Sebastian J. Vollmer, Piecewise Deterministic Markov Processes for Scalable Monte Carlo on Restricted Domains, arXiv preprint, 2017.\nJoris Bierkens, Paul Fearnhead and Gareth Roberts, The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data, arXiv preprint, 2016."
},

{
    "location": "examples/ex_gbps1.html#",
    "page": "Global BPS (Truncated Gaussian)",
    "title": "Global BPS (Truncated Gaussian)",
    "category": "page",
    "text": "(the code for this example can be found here)"
},

{
    "location": "examples/ex_gbps1.html#Global-BPS-(Truncated-Gaussian)-1",
    "page": "Global BPS (Truncated Gaussian)",
    "title": "Global BPS (Truncated Gaussian)",
    "category": "section",
    "text": "In this example we XXXXXXXStart by loading the library:using PDMPyou will then need to define two elements:a geometry (boundaries)\nan energy (gradient of the log-likelihood of the target)At the moment, the package can handle unconstrained geometries and polygonal domains (see XXXXX). Let's say we want to be constrained to the positive orthan in 2D:p = 2\n# normal to faces and intercepts\nns, a = eye(p), zeros(p)\ngeom  = Polygonal(ns, a)\n# for a given ray, which boundary does it hit?\nnextbd(x, v) = nextboundary(geom, x, v)Here ns and a are the normals and the intercepts of the faces. The type Polygonal encapsulates the geometry. The function nextboundary returns the next boundary on the current ray [x,x+tv] with t>0 as well as the time of the hit.We then need to specify a model: we need to define a function of the form gradll(x) which can return the gradient of the log-likelihood at some point x. Here, let us consider a 2D gaussian.# build a valid precision matrix, the cholesky decomposition of\n# the covariance matrix will be useful later to build a sensible\n# starting point.\nsrand(12)\nP1  = randn(p,p)\nP1 *= P1'\nP1 += norm(P1)/100*eye(p)\nC1  = inv(P1); C1 += C1'; C1/=2;\nL1  = cholfact(C1)\nmu  = zeros(p)+1.\nmvg = MvGaussianCanon(mu, P1)Here, we have defined the gaussian through the \"Canonical\" representation (see XXXXX) i.e.: by specifying a mean and a precision matrix.The gradient of the log-likelihood is then given bygradll(x) = gradloglik(mvg, x)Remark: if you want to implement your own model, you should define your model in (XXXXXX) and make sure it implements a gradloglik function.Next, we need to define the function which can return the first arrival time of the Inhomogenous Poisson Process (cf. algorithm). Note that you could be using nextevent_zz here as well if you wanted to use the Zig-Zag sampler (and you could implement other kernels as well, see HERE XXXXX).nextev(x, v) = nextevent_bps(mvg, x, v)For a Gaussian (and some other simple distributions), this is analytical through an inversion-like method (cf. algorithm). Another approach is the thinning approach using a bounding intensity. At the moment thinning with a linear bound is implemented (cf XXXXX).Finally, you need to specify the parameters of the simulation such as the starting point and velocity, the length of the path generated, the rate of refreshment and the maximum number of gradient evaluations. (see discussion)T    = 1000.0   # length of path generated\nlref = 2.0      # rate of refreshment\nx0   = mu+L1[:L]*randn(p) # sensible starting point\nv0   = randn(p) # starting velocity\nv0  /= norm(v0) # put it on the sphere (not necessary)\n# Define a simulation\nsim = Simulation( x0, v0, T, nextev, gradll,\n            nextbd, lref ; maxgradeval = 10000)And finally, generate the path and recover some details about the simulation.(path, details) = simulate(sim)The path object belongs to the type Path and can be sampled using samplepath.A crude test is to check that the estimated mean obtained through quadrature along the path yields a similar result as a basic Monte Carlo estimator.# Building a basic MC estimator\nsN = 1000\ns  = repmat(mu,1,sN)+L1[:L]*randn(p,sN)\nmt = zeros(2)\nnp = 0\n# Sum for all samples in the positive orthan\nss = [s; ones(sN)']\nmt = sum(ss[:,i] for i in 1:sN if !any(e->e<0, ss[1:p,i]))\nmt = mt[1:p]/mt[end]You can now compare the norm of mt to pathmean(path) and you will see that the relative error is below 5%."
},

{
    "location": "examples/ex_lbps1.html#",
    "page": "Local BPS (Chain of Gaussians)",
    "title": "Local BPS (Chain of Gaussians)",
    "category": "page",
    "text": "(the code for this example can be found here)"
},

{
    "location": "examples/ex_lbps1.html#Local-BPS-(Chain-of-Gaussians)-1",
    "page": "Local BPS (Chain of Gaussians)",
    "title": "Local BPS (Chain of Gaussians)",
    "category": "section",
    "text": "The approach to using the local BPS is much the same as for the global one except that you need to specify a FactorGraph. That object will contain the structure of the factor graph (which factor is connected to which variables) as well as the list of all factors (which have a gll and nextevent since each factor can be seen individually as a small BPS).Let's declare a chain of bivariate gaussians:using PDMP\nnfac = 3 # number of factors\n\nmvg             = PDMP.MvGaussianStandard(zeros(2),eye(2))\ngll(x)          = PDMP.gradloglik(mvg, x)\nnextevent(x, v) = PDMP.nextevent_bps(mvg, x, v)\n\n# all factors have that same likelihood\nchainfactor(i) = Factor(nextevent,gll,i)\n\n# assemble into a chain graph\nchain = chaingraph([chainfactor(i) for i in 1:nfac])This is a simple graph with a known structure so that it's already defined through the chaingraph function (in src/local/factorgraph.jl). For an arbitrary graph, you would need to provide two things:the structure of the factor graph: a list of list where each elementcorresponds to a factor and the corresponding list contains the indices of the variables attached to that factorthe list of factorsThe rest is very similar to the global BPS:lambdaref  = .01\nmaxnevents = 10000\nT          = Inf\nnvars      = chain.structure.nvars\nx0         = randn(nvars)\nv0         = randn(nvars)\nv0        /= norm(v0)\n\nlsim = LocalSimulation(chain, x0, v0, T, maxnevents, lambdaref)\n\n(all_evlist, details) = simulate(lsim)The all_evlist object contains a list of EventList corresponding to what happened on each of the factors. It can also be sampled using samplelocalpath (cf. src/local/event.jl)."
},

{
    "location": "techdoc/types.html#",
    "page": "Types",
    "title": "Types",
    "category": "page",
    "text": ""
},

{
    "location": "techdoc/types.html#PDMP.FactorGraph",
    "page": "Types",
    "title": "PDMP.FactorGraph",
    "category": "Type",
    "text": "FactorGraph\n\nEncapsulation of a factor graph. It is made out of a structure (FactorGraphStruct) and an array of factors corresponding to the given structure. See also: Factor, FactorGraphStruct.\n\n\n\n"
},

{
    "location": "techdoc/types.html#PDMP-Types-1",
    "page": "Types",
    "title": "PDMP Types",
    "category": "section",
    "text": "(/!\\WIP/!\\)FactorGraphInt for Int64 – in PDMP.jl\nFloat for Float64 – in PDMP.jl\nAllowedVarType for Union{Float, Vector{Float}} (variable types in the local BPS) – in local/event.jl\nAllowedTimeType for Union{Vector{Float}, LinSpace{Float}} (format of collections of times that can be passed to extract samples from events) – in path.jl"
},

{
    "location": "techdoc/types.html#Abstract-types-1",
    "page": "Types",
    "title": "Abstract types",
    "category": "section",
    "text": "MvGaussian for multivariate gaussians – in models/mvgaussian.jl\nDomain for domain description (in global BPS) – in geometry.jl\nIPPSamplingMethod for sampling methods of an IPP – in ippsampler.jl\nThinning <: IPPSamplingMethod for thinning methods of an IPP – in ippsampler.jl"
},

{
    "location": "techdoc/types.html#Specific-types-(global)-1",
    "page": "Types",
    "title": "Specific types (global)",
    "category": "section",
    "text": "geometryUnconstrained <: Domain (immutable) defines a domain without boundary (signature only)\nPolygonal <: Domain (immutable) defines a domain with affine boundaries determined by normals and interceptsippsamplerLinearBound <: Thinning (immutable) thinning method using a linear bound following a uniform bound on the eigenvalues of the Hessian\nNextEvent object returned when sampling from the IPP (contains a bouncing time, whether to bounce or not and a flipindex (ZZ case))pathPath container for a path of the global sampling (stores list of corners and times)\nSegment (immutable) container to store the two corners of a linear segment, makes it easier to sample from a PathsimulateSimulation container for the parameters of a global simulation"
},

{
    "location": "techdoc/types.html#Model-types-(global)-1",
    "page": "Types",
    "title": "Model types (global)",
    "category": "section",
    "text": "LogReg (immutable) model for a logistic regression"
},

{
    "location": "techdoc/types.html#Specific-types-(local)-1",
    "page": "Types",
    "title": "Specific types (local)",
    "category": "section",
    "text": "eventEvent{T<:AllowedVarType} (immutable) a triple (x,v,t) corresponding to an event for one of the node.\nEventList{T<:AllowedVarType} list of events (stored as lists of xs, vs, ts in order to be traversed more effectively than a Vector{Event}).\nAllEventList container for the EventList of all the nodesfactorgraphFactor (immutable) attaches the nextevent function (sampling from IPP), the gradient of the corresponding log-likelihood (gll) and an index\nFactorGraphStruct (immutable) contains the pattern of connections between nodes via flist and vlist (storing what variable is attached to what factor and vice-versa) nfactors and nvars keep track of the number of factors and the number of nodes (variables)\nFactorGraph (immutable) container for all the Factor and the FactorGraphStructsimulateLocalSimulation (immutable) container for the parameters of a local simulation (the factor graph, initial positions, etc)"
},

]}
