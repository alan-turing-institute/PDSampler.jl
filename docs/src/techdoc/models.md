# [Local sampler](@id td-models)

Link to the source files:

* [`models/mvgaussian.jl`](https://github.com/alan-turing-institute/PDSampler.jl/blob/master/src/models/mvgaussian.jl)
* [`models/logreg.jl`](https://github.com/alan-turing-institute/PDSampler.jl/blob/master/src/models/logreg.jl)
* [`models/pmf.jl`](https://github.com/alan-turing-institute/PDSampler.jl/blob/master/src/models/pmf.jl)

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

In the sequel we write $\mu$ the mean, $\Sigma$ the covariance matrix and $\Omega$ the precision matrix. The different way to parametrise the distributions are as follows:

* `MvGaussianStandard`, direct: $(\mu, \Sigma)$, indirect: (\Omega\mu,\Omega)
* `MvDiagonalGaussian`, direct: $(\mu, (\sigma_i))$, indirect: $(\sigma_i^2)$
* `MvGaussianCanon`, direct: $(\mu, \Omega)$, indirect: $(\Omega\mu)$
* `MvGaussianNatural`, direct: $(\Omega\mu,-\Omega)$

The preferred way is the "canonical" representation (most efficient).

**Note**: "direct" means that these are the parameters passed to the constructor while "indirect" means that these values are computed when the constructor is called.

### Auxiliary functions

Internally, the types mentioned above are shortened to `MvGS`, `MvDG` etc. Then a number of simplifying functions are defined (these simplify the computation of the log-likelihood and gradient of the log-likelihood)

* `mvg_mu` to recover $\mu$
* `mvg_precmu` to recover $\Omega\mu$
* `mvg_precmult` taking a point and multiplying it by $\Omega$

`gradloglik` is then trivial to compute.

## Logistic Regression

The logistic regression considers a feature matrix `X`, a response `y`, the Lipschitz constant associated to it and dimensionality parameters.

### Auxiliary functions

A number of auxiliary functions are defined to prevent numerical instabilities and ensure that the computation of the log-likelihood and gradient of the log-likelihood can be expressed simply.

The `gradloglik_cv` considers a control-variate gradient developed around a given point (see [this paper](https://arxiv.org/pdf/1701.04244.pdf) for more details).

**Note**: the response is in $\{-1,1\}$.

## Probabilistic Matrix Factorisation

This model considers a normal distribution on every entry of a matrix $r_{ij}$:

\begin{equation}
\mathcal N(r_{ij}; \langle u,v\rangle , \sigma^2)
\end{equation}

The resulting intensity can be shown to be a truncated cubic for which we can in fact also do exact sampling.

The `pmf_case*` correspond to the various possible cases depending on where the roots of the cubic are.
