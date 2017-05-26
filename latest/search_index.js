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
    "text": "PDMP.jl is a package designed to provide an efficient, flexible, and expandable framework for samplers based on Piecewise Deterministic Markov Processes and their applications. This includes the Bouncy Particle Sampler and the Zig-Zag Sampler. See the references at the bottom of this page.Pages = [\n    \"aboutpdmp.md\",\n    ]\nDepth = 1The package is currently under construction (Spring 2017). The project is hosted and maintained by the Alan Turing Institute (ATI) and currently developed by Thibaut Lienart. If you encounter problems, please open an issue on Github. If you have comments or wish to collaborate, please send an email to tlienart σ turing > ac > uk."
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
    "text": "The following examples will introduce you to the functionalities of the package.Pages = [\n    \"examples/ex_gbps1.md\",\n    \"examples/ex_lbps1.md\"\n    ]\nDepth = 1"
},

{
    "location": "index.html#Code-documentation-1",
    "page": "Introduction",
    "title": "Code documentation",
    "category": "section",
    "text": "These pages introduce you to the core of the package and its interface. This is useful if you are looking into expanding the code yourself to add a capacity or a specific model.Pages = [\n    \"techdoc/structure.md\",\n    \"techdoc/types.md\"\n    ]\nDepth = 1"
},

{
    "location": "index.html#Contributing-1",
    "page": "Introduction",
    "title": "Contributing",
    "category": "section",
    "text": "Pages = [\n    \"contributing/addingexample.md\",\n]\nDepth = 1"
},

{
    "location": "index.html#References-1",
    "page": "Introduction",
    "title": "References",
    "category": "section",
    "text": "Alexandre Bouchard-Côté, Sebastian J. Vollmer and Arnaud Doucet, The Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain Monte Carlo Method, arXiv preprint, 2015.\nJoris Bierkens, Alexandre Bouchard-Côté, Arnaud Doucet, Andrew B. Duncan, Paul Fearnhead, Gareth Roberts and Sebastian J. Vollmer, Piecewise Deterministic Markov Processes for Scalable Monte Carlo on Restricted Domains, arXiv preprint, 2017.\nJoris Bierkens, Paul Fearnhead and Gareth Roberts, The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data, arXiv preprint, 2016."
},

{
    "location": "aboutpdmp.html#",
    "page": "About PDMP",
    "title": "About PDMP",
    "category": "page",
    "text": ""
},

{
    "location": "aboutpdmp.html#About-PDMP-samplers-1",
    "page": "About PDMP",
    "title": "About PDMP samplers",
    "category": "section",
    "text": "This page aims at giving a very brief introduction to the concept of PDMP samplers (below we will refer to the algorithm but it should be understood as a class of algorithms). We also give some insight into how it is implemented although we cover the implementation in more details in the technical documentation. This is not meant to be a rigorous presentation of the algorithm (for this, please see the references at the bottom of this page). Rather, we focus here on the \"large building blocks\" behind the algorithm."
},

{
    "location": "aboutpdmp.html#Basic-idea-(global-samplers)-1",
    "page": "About PDMP",
    "title": "Basic idea (global samplers)",
    "category": "section",
    "text": "The purpose of the algorithm is to be able to evaluate expected values with respect to an arbitrary target distribution which we assume admits a probability density function  pi. For simplicity, we assume that piCto mathbb R^+ with Csubseteq mathbb R^p, convex. The objective is therefore to compute a weighted integral of the form:\\begin{equation}     \\mathbb E_{\\pi}[\\varphi(X)] = \\int_{C} \\varphi(x)\\pi(x)\\,\\mathrm{d}x \\end{equation}The samples considered generate a piecewise-linear path\\begin{equation}     x(t) = x^{(i)} + v^{(i)}(t-t_i) \\quad \\text{for}\\quad t\\in[t_i, t_{i+1}] \\end{equation}determined by an initial position x^(0) and velocity v^(0) at time t_0=0 and a set of positive event times t_1t_2dots. Under some conditions for the generation of the times and the velocities, the expected value can be approximated with\\begin{eqnarray}     \\mathbb E_{\\pi}[\\varphi(X)] &\\approx& {1\\over T} \\int_0^T\\varphi(x(t))\\mathrm{d}t \\end{eqnarray}and the integral in the right hand side can be expressed as a sum of one-dimensional integrals along each linear segment."
},

{
    "location": "aboutpdmp.html#Generating-times-and-velocities-1",
    "page": "About PDMP",
    "title": "Generating times and velocities",
    "category": "section",
    "text": "The algorithm generates a sequence of triples of the form (t^(i) x^(i) v^(i)). Let us assume that the algorithm is currently at one of those event points and show how to compute the next triple. To do so, the algorithm executes the following steps:it generates a travel time tau drawing from a specific process,\nthe next position is obtained by traveling along the current ray for the travel time tau i.e.: x^(i+1) = x^(i) + tau v^(i),\na new velocity v^(i+1) is generated.First we will explore how the travel time is generated and then how the new velocity is computed."
},

{
    "location": "aboutpdmp.html#Sampling-a-travel-time-1",
    "page": "About PDMP",
    "title": "Sampling a travel time",
    "category": "section",
    "text": "The travel time tau is obtained as the minimum of three times which we will denote by tau_b tau_h tau_r. Following the case, the computation of the new velocity will be different.The first (and most important) one, tau_b, is the first arrival time of an Inhomogenous Poisson Process (IPP) with an intensity that should verify some properties with respect to the target distribution. The Bouncy Particle Sampler (BPS) in particular considers the following intensity with U the log-likelihood of the (unnormalised) target pi:\\begin{eqnarray}     \\lambda(\\tau; x, v) = \\langle \\nabla U(x + \\tau v ), v \\rangle^+ \\end{eqnarray}where x and v are the current points and f^+=max(f0). Sampling from an IPP is not trivial in general but there are a few well known techniques that can be applied depending on the target.The other two times are easy to compute:the first, tau_h, is the time of first hit with the boundary of the domain C along the current ray x(t)=x^(i)+(t-t_i)v^(i) for tt_i. This guarantees that the trajectory stays in C.\nthe second, tau_r, is a refreshment time sampled from an exponential distribution with a fixed rate. This guarantees full exploration of C (see BPS paper for details)."
},

{
    "location": "aboutpdmp.html#Computing-a-new-velocity-(BPS)-1",
    "page": "About PDMP",
    "title": "Computing a new velocity (BPS)",
    "category": "section",
    "text": "Below we discuss the case of the BPS, the computations can be different for different samplers (such as the ZZ) but the essence of the method is the same.As mentioned above, we take tau = min(tau_b tau_h tau_r). Depending on the case, three actions can be takena bounce with tau = tau_b where the new velocity is obtained by specular reflection against the tangent to the gradient of the log-likelihood at the point x(tau_b),\na boundary bounce with tau=tau_h where the new velocity is obtained by specular reflection against the tangent to the boundary at the point of hit x(tau_h),\na refreshment with tau=tau_r where the new velocity is drawn from a reference process such as a spherical Gaussian.The update of the velocity goes as follows for the BPS (specular reflection):\\begin{equation}     v \\leftarrow v - 2\\langle \\nabla U(x), v\\rangle{\\nabla U(x)\\over \\|\\nabla U(x)\\|^2}. \\end{equation}The figure below illustrates the specular reflexion, starting at the red point and going along the current ray (red, dashed line), we have a new event corresponding to a bounce or a hit (blue dot). In both case a specular reflection is executed (blue dashed line). The black line represents the tangent to either the boundary at that point or to the log-likelihood depending on the case.(Image: )"
},

{
    "location": "aboutpdmp.html#Putting-the-pieces-together-1",
    "page": "About PDMP",
    "title": "Putting the pieces together",
    "category": "section",
    "text": "The simple global sampler can be expressed as follows:Initialize (x^(0) v^(0)) and T the trajectory length\nFor i=12dots, consider the ray x^(i-1)+tau v^(i-1)\nSimulate tau_b from an IPP along the ray\nCompute tau_h, simulate tau_r and let tau=min(tau_htau_rtau_b)\nFollowing the case in (2.) compute the new velocity v^(i)\nStore the new triple (t_i-1+tau x^(i-1)+tau v^(i-1) v^(i))\nif t_i ge T stop.\nReturn the path: (t_i x^(i) v^(i))_i=01dotsFollowing this representation, here are the key files of the code:A way to sample from an IPP:  ippsampler.jl.\nA way to define the geometry and in particular to compute the next boundary hit when traveling along a given ray: geometry.jl.\nA way to define how the velocity needs to be updated (reflection, refreshments): kernels.jl.\nA way to store a path formed of triples and compute integrals along it: path.jl.\nA core loop: simulate.jl.We describe those in details and give explanations as to how to expand the toolbox in the technical documentation part."
},

{
    "location": "aboutpdmp.html#Local-Samplers-1",
    "page": "About PDMP",
    "title": "Local Samplers",
    "category": "section",
    "text": ""
},

{
    "location": "aboutpdmp.html#References-1",
    "page": "About PDMP",
    "title": "References",
    "category": "section",
    "text": "Alexandre Bouchard-Côté, Sebastian J. Vollmer and Arnaud Doucet, The Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain Monte Carlo Method, arXiv preprint, 2015.\nJoris Bierkens, Alexandre Bouchard-Côté, Arnaud Doucet, Andrew B. Duncan, Paul Fearnhead, Gareth Roberts and Sebastian J. Vollmer, Piecewise Deterministic Markov Processes for Scalable Monte Carlo on Restricted Domains, arXiv preprint, 2017.\nJoris Bierkens, Paul Fearnhead and Gareth Roberts, The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data, arXiv preprint, 2016."
},

{
    "location": "examples/ex_gbps1.html#",
    "page": "Global BPS",
    "title": "Global BPS",
    "category": "page",
    "text": ""
},

{
    "location": "examples/ex_gbps1.html#Global-BPS-(Truncated-Gaussian)-1",
    "page": "Global BPS",
    "title": "Global BPS (Truncated Gaussian)",
    "category": "section",
    "text": "(the code for this example can be found here, note that the doc rendered here was automatically generated, if you want to fix it, please do it in the julia code directly*)In this example we XXXXXXXStart by loading the library:using PDMPyou will then need to define two elements:a geometry (boundaries)\nan energy (gradient of the log-likelihood of the target)At the moment, the package can handle unconstrained geometries and polygonal domains (see XXXXX). Let's say we want to be constrained to the positive orthan in 2D:p = 2\n# normal to faces and intercepts\nns, a = eye(p), zeros(p)\ngeom  = Polygonal(ns, a)\n# for a given ray, which boundary does it hit?\nnextbd(x, v) = nextboundary(geom, x, v)Here ns and a are the normals and the intercepts of the faces. The type Polygonal encapsulates the geometry. The function nextboundary returns the next boundary on the current ray [x,x+tv] with t>0 as well as the time of the hit.We then need to specify a model: we need to define a function of the form gradll(x) which can return the gradient of the log-likelihood at some point x. Here, let us consider a 2D gaussian.# build a valid precision matrix, the cholesky decomposition of\n# the covariance matrix will be useful later to build a sensible\n# starting point.\nsrand(12)\nP1  = randn(p,p)\nP1 *= P1'\nP1 += norm(P1)/100*eye(p)\nC1  = inv(P1); C1 += C1'; C1/=2;\nL1  = cholfact(C1)\nmu  = zeros(p)+1.\nmvg = MvGaussianCanon(mu, P1)Here, we have defined the gaussian through the \"Canonical\" representation (see XXXXX) i.e.: by specifying a mean and a precision matrix.The gradient of the log-likelihood is then given bygradll(x) = gradloglik(mvg, x)Remark: if you want to implement your own model, you should define your model in (XXXXXX) and make sure it implements a gradloglik function.Next, we need to define the function which can return the first arrival time of the Inhomogenous Poisson Process (cf. algorithm). Note that you could be using nextevent_zz here as well if you wanted to use the Zig-Zag sampler (and you could implement other kernels as well, see HERE XXXXX).nextev(x, v) = nextevent_bps(mvg, x, v)For a Gaussian (and some other simple distributions), this is analytical through an inversion-like method (cf. algorithm). Another approach is the thinning approach using a bounding intensity. At the moment thinning with a linear bound is implemented (cf XXXXX).Finally, you need to specify the parameters of the simulation such as the starting point and velocity, the length of the path generated, the rate of refreshment and the maximum number of gradient evaluations. (see discussion)T    = 1000.0   # length of path generated\nlref = 2.0      # rate of refreshment\nx0   = mu+L1[:L]*randn(p) # sensible starting point\nv0   = randn(p) # starting velocity\nv0  /= norm(v0) # put it on the sphere (not necessary)\n# Define a simulation\nsim = Simulation( x0, v0, T, nextev, gradll,\n            nextbd, lref ; maxgradeval = 10000)And finally, generate the path and recover some details about the simulation.(path, details) = simulate(sim)The path object belongs to the type Path and can be sampled using samplepath.A crude test is to check that the estimated mean obtained through quadrature along the path yields a similar result as a basic Monte Carlo estimator.# Building a basic MC estimator\nsN = 1000\ns  = repmat(mu,1,sN)+L1[:L]*randn(p,sN)\nmt = zeros(2)\nnp = 0\n# Sum for all samples in the positive orthan\nss = [s; ones(sN)']\nmt = sum(ss[:,i] for i in 1:sN if !any(e->e<0, ss[1:p,i]))\nmt = mt[1:p]/mt[end]You can now compare the norm of mt to pathmean(path) and you will see that the relative error is below 5%."
},

{
    "location": "examples/ex_lbps1.html#",
    "page": "Local BPS",
    "title": "Local BPS",
    "category": "page",
    "text": ""
},

{
    "location": "examples/ex_lbps1.html#Local-BPS-(Chain-of-Gaussians)-1",
    "page": "Local BPS",
    "title": "Local BPS (Chain of Gaussians)",
    "category": "section",
    "text": "(the code for this example can be found here, note that the doc rendered here was automatically generated, if you want to fix it, please do it in the julia code directly*)The approach to using the local BPS is much the same as for the global one except that you need to specify a FactorGraph. That object will contain the structure of the factor graph (which factor is connected to which variables) as well as the list of all factors (which have a lgradll and nextevent since each factor can be seen individually as a small BPS).Let's declare a chain of bivariate gaussians:using PDMP\nnfac = 3 # number of factors\n\nmvg = MvGaussianStandard(zeros(2),eye(2))\n\n# all factors have that same likelihood\nchainfactor(i) = Factor(\n                    (x,v)->nextevent_bps(mvg, x, v),\n                    (x)->gradloglik(mvg, x),\n                    i )\n\n# assemble into a chain graph\nchain = chaingraph([chainfactor(i) for i in 1:nfac])This is a simple graph with a known structure so that it's already defined through the chaingraph function (in src/local/factorgraph.jl). For an arbitrary graph, you would need to provide two things:the structure of the factor graph: a list of list where each elementcorresponds to a factor and the corresponding list contains the indices of the variables attached to that factorthe list of factorsThe rest is very similar to the global BPS:srand(123)\nlambdaref  = .01\nmaxnevents = 10000\nT          = Inf\nnvars      = chain.structure.nvars\nx0         = randn(nvars)\nv0         = randn(nvars)\nv0        /= norm(v0)\n\nlsim = LocalSimulation(chain, x0, v0, T, maxnevents, lambdaref)\n\n(all_evlist, details) = simulate(lsim)The all_evlist object contains a list of EventList corresponding to what happened on each of the factors. It can also be sampled using samplelocalpath (cf. src/local/event.jl)."
},

{
    "location": "techdoc/structure.html#",
    "page": "Code structure",
    "title": "Code structure",
    "category": "page",
    "text": ""
},

{
    "location": "techdoc/structure.html#Structure-of-the-code-1",
    "page": "Code structure",
    "title": "Structure of the code",
    "category": "section",
    "text": "blah blah"
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

{
    "location": "contributing/addingexample.html#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "contributing/addingexample.html#Adding-an-example-1",
    "page": "Examples",
    "title": "Adding an example",
    "category": "section",
    "text": "Examples are a great way to showcase how to make use of a specific feature. We consider two types of examples:a simple example (running in < 10 seconds)\na complex example (the complement of the first category)The first category is great as they can be used as tests and as documentation.The process is somewhat automated here and essentially all you have to do is write the example in the test directory and comment it accordingly, we show this below."
},

{
    "location": "contributing/addingexample.html#Syntax-for-the-example-1",
    "page": "Examples",
    "title": "Syntax for the example",
    "category": "section",
    "text": "Let's say you have an example which can be run in a few seconds and uses a new feature. You can effectively use it as a unit test by itself. To respect conventions, please name your example ex_NAME.jl and put it in test/.Start your code byusing Base.TestThen a few markers should be considered:#@startexample NAME_OF_EXAMPLE indicates that you start the code of the example proper, everything between that mark and the end flag will appear in the doc.\nEncapsulates all the explanations you want to appear in markdown between #= and =#, all the other comments will be taken as part of the julia code and shown in code blocks.\n#@endexample to indicate that the example is finished\nwrite a few tests that check that the example produces the right answer (unit test)So it should look like (a full example can be seen here)using Base.Test\n#@startexample A simple example\n#=\nIn this example we will show how to find the maximum over a collection of\nrandom numbers in Julia for all the numbers that are *negative*.\n=#\na = randn(500)\n# we use a comprehension\nm = maximum(a[i] for i in 1:length(a) if a[i]<0)\n#=\nThat's it!\n=#\n#@endexample\n@test m < 0.0Make sure the tests pass! This will generate the following markdown (see next point for the command that generates it):# Name of the Test\n\nIn this example we will show how to find the maximum over a collection of\nrandom numbers in Julia for all the numbers that are negative.\n ```julia\na = randn(500)\n# we use a comprehension\nm = maximum(a[i] for i in 1:length(a) if a[i]<0)\n ```\nThat's it!Remark, the spaces in front of the triple backticks in the markdown snippet above are not actually generated when you use readexamples.jl. The spaces are used here in order to escape these triple backticks so that the snippet does not end up being unduly fragmented in three pieces."
},

{
    "location": "contributing/addingexample.html#Declaring-your-example-1",
    "page": "Examples",
    "title": "Declaring your example",
    "category": "section",
    "text": "You have to mention your example in a few spots:in test/runtests.jl, add a line at the bottom following the examples already present i.e. something like @testset \"ex_NAME\"    begin include(\"ex_NAME.jl\") end (make sure this passes!)\nin docs/src/make.jl, add a line under \"Examples\" following the syntax of examples already present so something like: \"Name of your expl\" => \"examples/ex_name_of_example.jl\"Finally, to generate the markdown run the following command (this will act for all the examples at once so it will also refresh any other modification you may have added to other examples):julia docs/readexamples.jl"
},

]}
