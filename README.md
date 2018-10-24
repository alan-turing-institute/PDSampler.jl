# PDSampler

Unix | Windows | CodeCov | Docs | License
---- | ------- | ------- | ---- | -------
[![Travis](https://travis-ci.org/alan-turing-institute/PDSampler.jl.svg?branch=master)](https://travis-ci.org/alan-turing-institute/PDSampler.jl) | [![AppVeyor](https://ci.appveyor.com/api/projects/status/github/alan-turing-institute/PDSampler.jl?branch=master&svg=true)](https://ci.appveyor.com/project/tlienart/PDSampler-jl) | [![CodeCov](http://codecov.io/github/alan-turing-institute/PDSampler.jl/coverage.svg?branch=master)](http://codecov.io/github/alan-turing-institute/PDSampler.jl?branch=master) | [![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://alan-turing-institute.github.io/PDSampler.jl/latest) | [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

PDSampler.jl is a package designed to provide an efficient, flexible, and expandable framework for samplers based on *Piecewise Deterministic Markov Processes* and their applications.
This includes the **Bouncy Particle Sampler** and the **Zig-Zag Sampler**.

Please refer to [**the documentation**](https://alan-turing-institute.github.io/PDSampler.jl/latest) for information on how to use/expand this package.
The project is hosted by the [Alan Turing Institute](https://www.turing.ac.uk) (ATI). If you encounter problems, please [open an issue on Github](https://github.com/alan-turing-institute/PDSampler.jl/issues).
If you have comments or wish to collaborate, please send an email to `tlienart > cpg σ gmail > com`.

If you find this toolbox useful please star the repo. If you use it in your work, please cite this code and send us an email so that we can cite your work here.

If you want to make suggestions, if you want new features, please don't hesitate, open an issue or send an email.

### Contributors

* [Thibaut Lienart](http://www.stats.ox.ac.uk/~lienart/) (main dev)
* [Sebastian Vollmer](https://www2.warwick.ac.uk/fac/sci/maths/people/staff/vollmer/)
* [Andrew Duncan](http://www.imperial.ac.uk/complex-multiscale-systems/our-group/former-members/dr-andrew-duncan/)
* Martin O'Reilly (ATI)

## Installation and requirements

(This is explained in more details in the documentation)

Requirements:

* Julia ∈ `[0.7.*, 1.0.*]`, if you're on `0.6`, check out [the last legacy release](https://github.com/alan-turing-institute/PDSampler.jl/releases/tag/v0.1).

In the Julia REPL:

```julia
] add PDSampler
using PDSampler
```

Note that loading the package may take several seconds as some of the dependencies (in particular [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl) are quite slow to load).

## References

* Alexandre Bouchard-Côté, Sebastian J. Vollmer and Arnaud Doucet, [*The Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain Monte Carlo Method*](https://arxiv.org/abs/1510.02451), arXiv preprint, 2015.
* Joris Bierkens, Alexandre Bouchard-Côté, Arnaud Doucet, Andrew B. Duncan, Paul Fearnhead, Gareth Roberts and Sebastian J. Vollmer, [*Piecewise Deterministic Markov Processes for Scalable Monte Carlo on Restricted Domains*](https://arxiv.org/pdf/1701.04244.pdf), arXiv preprint, 2017.
* Joris Bierkens, Paul Fearnhead and Gareth Roberts, [*The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data*](https://arxiv.org/pdf/1607.03188.pdf), arXiv preprint, 2016.
* Changye Wu, Christian Robert, [*Generalized Bouncy Particle Sampler*](https://arxiv.org/pdf/1706.04781.pdf), arXiv preprint, 2017.

*Note*: if your paper is not listed here and you feel like it should, please open an issue (same goes if there is a mistake or if a preprint is now a proper-print).
