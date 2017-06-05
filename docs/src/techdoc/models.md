# [Local sampler](@id td-models)

Link to the source files:

* [`mvgaussian.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/models/mvgaussian.jl)
* [`logreg.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/models/logreg.jl)
* [`mvgaussian.jl`](https://github.com/alan-turing-institute/PDMP.jl/blob/master/src/models/pmf.jl)

## Generic model

A model must be an immutable type with an associated `gradloglik` function. It is important this function be coded *as efficiently as possible* since it is called a large number of time in any simulation.

## Multivariate Gaussian

### Hierarchy of types

Multiple parametrisation are possible. Some being more efficient than others while some are maybe more intuitive than others.

```
MvGaussian (abstract)
| — MvGaussianStandard
| — MvDiagonalGaussian
| — MvGaussianCanon
| — MvGaussianNatural
```

In the sequel we write $\mu$ the mean, $\Sigma$ the covariance matrix and $\Omega$ the precision matrix. The parametrisation are then as follows:

* `MvGaussianStandard`: direct: $(\mu, \Sigma)$, indirect: (\Omega\mu,\Omega)
* `MvDiagonalGaussian`: direct: $(\mu, \sigma)$

**Note**: "direct" means that these are the parameters passed to the constructor while "indirect" means that these values are computed when the constructor is called.

### Auxiliary functions



## Logistic Regression

## Probabilistic Matrix Factorisation
