# PDMP

Unix | Windows | CodeCov | Docs | License
---- | ------- | ------- | ---- | -------
[![Travis](https://travis-ci.org/alan-turing-institute/PDMP.jl.svg?branch=master)](https://travis-ci.org/alan-turing-institute/PDMP.jl) | [![AppVeyor](https://ci.appveyor.com/api/projects/status/github/alan-turing-institute/PDMP.jl?branch=master&svg=true)](https://ci.appveyor.com/project/tlienart/pdmp-jl) | [![CodeCov](http://codecov.io/github/alan-turing-institute/PDMP.jl/coverage.svg?branch=master)](http://codecov.io/github/alan-turing-institute/PDMP.jl?branch=master) | [![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://alan-turing-institute.github.io/PDMP.jl/latest) | [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

PDMP.jl is a package designed to provide an efficient, flexible, and expandable framework for samplers based on *Piecewise Deterministic Markov Processes* and their applications.
This includes the **Bouncy Particle Sampler** and the **Zig-Zag Sampler**.

Please refer to [**the documentation**](https://alan-turing-institute.github.io/PDMP.jl/latest) for information on how to use/expand this package. The project is hosted and maintained by the [Alan Turing Institute](https://www.turing.ac.uk) (ATI) and currently developed by [Thibaut Lienart](http://www.stats.ox.ac.uk/~lienart/). If you encounter problems, please [open an issue on Github](https://github.com/alan-turing-institute/PDMP.jl/issues).
If you have comments or wish to collaborate, please send an email to `tlienart σ turing > ac > uk`.

If you find this toolbox useful please star the repo. If you use it in your work, please cite this code and send us an email so that we can cite your work here.

## Installation and requirements

(This is explained in more details in the documentation)

Requirements:

* Julia in `[0.5.x, 0.6.x]` (*Julia `0.7` introduces breaking changes to the syntax, since it's still a "nightly" release we are sticking to rc versions.*)
* 64-bit architecture (*It may be possible to adapt the code to run on 32 bit but it would be a bit of a headache, if there is genuine interest for this, please open issue and explain.*)

In the Julia REPL:

```julia
Pkg.clone("https://github.com/alan-turing-institute/PDMP.jl.git")
using PDMP
```

## References

* Alexandre Bouchard-Côté, Sebastian J. Vollmer and Arnaud Doucet, [*The Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain Monte Carlo Method*](https://arxiv.org/abs/1510.02451), arXiv preprint, 2015.
* Joris Bierkens, Alexandre Bouchard-Côté, Arnaud Doucet, Andrew B. Duncan, Paul Fearnhead, Gareth Roberts and Sebastian J. Vollmer, [*Piecewise Deterministic Markov Processes for Scalable Monte Carlo on Restricted Domains*](https://arxiv.org/pdf/1701.04244.pdf), arXiv preprint, 2017.
* Joris Bierkens, Paul Fearnhead and Gareth Roberts, [*The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data*](https://arxiv.org/pdf/1607.03188.pdf), arXiv preprint, 2016.
* Changye Wu, Christian Robert, [*Generalized Bouncy Particle Sampler*](https://arxiv.org/pdf/1706.04781.pdf), arXiv preprint, 2017.

*Note*: if your paper is not listed here and you feel like it should, please open an issue (same goes if there is a mistake or if a preprint is now a proper-print).
