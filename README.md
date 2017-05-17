# PDMP

Unix | Windows | CodeCov | Docs | License
---- | ------- | ------- | ---- | -------
[![Travis](https://travis-ci.org/alan-turing-institute/PDMP.jl.svg?branch=master)](https://travis-ci.org/alan-turing-institute/PDMP.jl) | [![AppVeyor](https://ci.appveyor.com/api/projects/status/github/alan-turing-institute/PDMP.jl?branch=master&svg=true)](https://ci.appveyor.com/project/tlienart/pdmp-jl) | [![CodeCov](http://codecov.io/github/alan-turing-institute/PDMP.jl/coverage.svg?branch=master)](http://codecov.io/github/alan-turing-institute/PDMP.jl?branch=master) | [![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://alan-turing-institute.github.io/PDMP.jl/latest) | [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

PDMP.jl is a package designed to provide an efficient, flexible, and expandable framework for samplers based on *Piecewise Deterministic Markov Processes* and their applications.
This includes the **Bouncy Particle Sampler** and the **Zig-Zag Sampler**.

Please refer to [**the documentation**](https://alan-turing-institute.github.io/PDMP.jl/latest) for information on how to use/expand this package.

## Installation and requirements

(This is explained in more details in the documentation)

Requirements:

* Julia >= 0.5
* 64-bit architecture

In the Julia REPL:

```julia
Pkg.clone("https://github.com/alan-turing-institute/PDMP.jl.git")
using PDMP
```
