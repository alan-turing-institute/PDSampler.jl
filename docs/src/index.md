# PDMP.jl Documentation

[PDMP.jl](https://github.com/alan-turing-institute/PDMP.jl) is a package designed to provide an efficient, flexible, and expandable framework for samplers based on *Piecewise Deterministic Markov Processes* and their applications.
This includes the **Bouncy Particle Sampler** and the **Zig-Zag Sampler**.
See [the references](#references) at the bottom of this page.

```@contents
Pages = [
    "aboutpdmp.md",
    ]
Depth = 2
```

The package is currently under construction (Spring 2017).
The project is hosted and maintained by the [Alan Turing Institute](https://www.turing.ac.uk) (ATI) and currently developed by [Thibaut Lienart](http://www.stats.ox.ac.uk/~lienart/). If you encounter problems, please [open an issue on Github](https://github.com/alan-turing-institute/PDMP.jl/issues).
If you have comments or wish to collaborate, please send an email to `tlienart σ turing > ac > uk`.

### Using the Package

To install the (currently unregistered) package, use the following command inside the Julia REPL:

```julia
Pkg.clone("git://github.com/alan-turing-institute/PDMP.jl.git")
```

To load the package, use the command:

```julia
using PDMP
```

You can also run the tests with ```Pkg.test("PDMP")``` and update to the latest Github version with ```Pkg.update("PDMP")```.

### Examples

The following examples will introduce you to the functionalities of the package.

```@contents
Pages = [
    "examples/ex_gbps1.md",
    "examples/ex_lbps1.md"
    ]
Depth = 2
```

### Code documentation

These pages introduce you to the core of the package and its interface.
This is useful if you are looking into expanding the code yourself to add a capacity or a specific model.

```@contents
Pages = [
    "techdoc/structure.md",
    "techdoc/types.md"
    ]
Depth = 2
```

### Contributing

```@contents
Pages = [
    "contributing/addingexample.md",
]
Depth = 2
```

### References

* Alexandre Bouchard-Côté, Sebastian J. Vollmer and Arnaud Doucet, [*The Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain Monte Carlo Method*](https://arxiv.org/abs/1510.02451), arXiv preprint, 2015.
* Joris Bierkens, Alexandre Bouchard-Côté, Arnaud Doucet, Andrew B. Duncan, Paul Fearnhead, Gareth Roberts and Sebastian J. Vollmer, [*Piecewise Deterministic Markov Processes for Scalable Monte Carlo on Restricted Domains*](https://arxiv.org/pdf/1701.04244.pdf), arXiv preprint, 2017.
* Joris Bierkens, Paul Fearnhead and Gareth Roberts, [*The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data*](https://arxiv.org/pdf/1607.03188.pdf), arXiv preprint, 2016.
