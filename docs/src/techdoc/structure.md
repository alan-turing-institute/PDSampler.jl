# Code structure

In this part, we discuss briefly how the code is organised and the role of the key files as well as the workflow for extending the code.

## General notes

A few design choices have been made and should be respected (or modified with a good reason and this section scrapped):

* `Float` stands for `Float64` assuming that everything is done on 64-bit architecture.
* When possible, abstract types are created to suggest a hierarchy of types. This helps readability and generalisation (see for example in `ippsampler.jl`, abstract type: `IPPSamplingMethod` and `Thinning <: IPPSamplingMethod` and `LinearBound <: Thinning`)

## Source files

The structure of the `src/` folder is as follows:

```
├── PDMP.jl
├── geometry.jl
├── ippsampler.jl
├── kernels.jl
├── local
│   ├── event.jl
│   ├── factorgraph.jl
│   └── simulate.jl
├── models
│   ├── logreg.jl
│   ├── mvgaussian.jl
│   └── pmf.jl
├── path.jl
└── simulate.jl
```

The central file is `PDMP.jl` which serves one key purpose: declaring what the package needs (`Compat`, `Polynomials`, ...) and including the files that contain the effective pieces of code. It also *exports* some generic functions that are used throughout the package.

**Note**: in Julia everything should be wrapped around by a *module*. The `using PkgName` indicates that we want to have access to the functions exported by the package `PkgName` in the current scope (e.g.: the scope of the wrapping module or that of the REPL). The `export functionName` indicates that if another user wants to use our module (by entering `using PDMP`) s/he will have access to all of those functions directly.

Here is a high-level overview of the rest of the folder structure:

* `geometry`, `ippsampler`, `kernels` ([specific documentation](@ref td-coretools)): generic tools used throughout the package
* `path`, `simulate` ([specific documentation](@ref td-globalsampler)): tools to describe the path and how the simulation is run in the global case.
* `models/*` ([specific documentation](@ref td-models)): to define specific models, their likelihood, gradient of log-likelihood etc.
* `local/*` ([specific documentation](@ref td-localsampler)): to define events, factor graphs and how to run the algorithm in the local case.

## Test files

The `test/` folder contains a number of test files (one for each source file and one per executable example):

```
├── ex_gbps1.jl
├── ex_lbps1.jl
├── gaussian_test.jl
├── geometry_test.jl
├── ippsampler_test.jl
├── kernels_test.jl
├── local_event_test.jl
├── local_factorgraph_test.jl
├── local_simulate_test.jl
├── logreg_test.jl
├── path_test.jl
├── pmf_test.jl
├── runtests.jl
└── simulate_test.jl
```

Note that a few start with `ex_` these are *executable examples* which also serve as partial tests **and** as documentation.
The philosophy here is to have as many tests as possible that would break if anything is introduced in the code that could break other parts. These tests are not perfect and some may indeed need to be reinforced/fixed but at least provide some safeguards against harmful code modifications.

## Doc files

The `docs/` folder contains a large number of files. The part that is of interest is represented below:

```
├── build
│   ├── ...
├── make.jl
├── readexamples.jl
├── site
│   ├── ...
└── src
    ├── aboutpdmp.md
    ├── assets
    │   ├── ...
    ├── contributing
    │   ├── addingexample.md
    │   └── addingfeature.md
    ├── examples
    │   ├── ex_gbps1.md
    │   └── ex_lbps1.md
    ├── index.md
    └── techdoc
        ├── coretools.md
        ├── global.md
        ├── local.md
        ├── models.md
        ├── structure.md
        └── types.md
```

The `make.jl` file is the central file which dictates how the documentation is to be built. It can be executed in a Julia REPL (provided you have added the `Documenter` package) and you can then locally see the updated version of the documentation by opening `build/index.html`.
The `readexamples.jl` file transforms the example files `test/ex_*` into publishable examples.

**Note**: if you are editing the documentation and wish to compile it, the recommendation is to keep your REPL open. The first compiling will be a bit slow (Documenter warming up) the next ones will be pretty much instantaneous with possibly a lot of warning messages about docstrings not having been found for every function, you can safely ignore all of that and just refresh the page `build/index.html` in your browser.
