export
    MvGaussian,
    MvGaussianStandard,
    MvDiagonalGaussian,
    MvGaussianCanon,
    MvGaussianNatural

abstract type MvGaussian end

"""
    MvGaussianStandard(mu, cov)

Standard representation of a multivariate Gaussian instantiated from the mean
and the covariance matrix. Users should prefer the Canonical representation if
they can use that.
"""
struct MvGaussianStandard <: MvGaussian
    mu::Vector{AbstractFloat}  # Precision*mean
    cov::Matrix{AbstractFloat} # negative Precision
    p::Int # Dimensions
    #
    prec::Matrix{AbstractFloat}
    precmu::Vector{AbstractFloat}
    function MvGaussianStandard(mu, cov)
        # stored+considered as lower triangular
        L = cholfact(Symmetric(cov))[:L]
        Li = inv(L)
        # going through CholFact is more stable + guarantees Positive Def
        prec = Li' * Li
        new(mu, cov, length(mu), prec, prec * mu)
    end
end

"""
    MvDiagonalGaussian(mu,sigma)

Representation of a multivariate diagonal gaussian (diagonal covariance).
"""
struct MvDiagonalGaussian <: MvGaussian
    mu::Vector{AbstractFloat}     # mean
    sigma::Vector{AbstractFloat}  # standard deviation
    sigma2::Vector{AbstractFloat} # standard deviation squared
    p::Int # dimensions
    MvDiagonalGaussian(mu, sigma) = new(mu, sigma, sigma.^2, length(mu))
end

"""
    MvGaussianCanon

Canonical representation of a multivariate Gaussian instantiated from the mean
and the precision matrix. This is the preferred representation.
"""
struct MvGaussianCanon <: MvGaussian
    mu::Vector{AbstractFloat}     # mean
    prec::Matrix{AbstractFloat}   # Precision
    precmu::Vector{AbstractFloat} # Precision*mean
    p::Int # Dimensions
    function MvGaussianCanon(mu, prec)
        @assert length(mu) == size(prec, 1) == size(prec, 2)
        new(mu, prec, prec * mu, length(mu))
    end
end

"""
    MvGaussianNatural

Natural parameter space representation of a multivariate Gaussian.
"""
struct MvGaussianNatural <: MvGaussian
    precmu::Vector{AbstractFloat}   # Precision*mean
    negprec::Matrix{AbstractFloat}  # negative Precision
    p::Int # Dimensions
    MvGaussianNatural(precmu, negprec) = new(precmu, negprec, length(precmu))
end

# -----------------------------------------------------------------------------
### Helper functions for efficient computations
# the type aliases are not exposed

const MvGS = MvGaussianStandard
const MvGC = MvGaussianCanon
const MvGN = MvGaussianNatural
const MvDG = MvDiagonalGaussian

mvg_mu(g::MvGC) = g.mu
mvg_mu(g::MvGN) = -g.negprec \ g.precmu
mvg_mu(g::MvGS) = g.mu
mvg_mu(g::MvDG) = g.mu

mvg_precmu(g::MvGC) = g.precmu
mvg_precmu(g::MvGN) = g.precmu
mvg_precmu(g::MvGS) = g.precmu
mvg_precmu(g::MvDG) = g.mu ./ g.sigma2

mvg_precmult(g::MvGC, x::Vector{AbstractFloat}) =  g.prec * x
mvg_precmult(g::MvGN, x::Vector{AbstractFloat}) = -g.negprec * x
mvg_precmult(g::MvGS, x::Vector{AbstractFloat}) =  g.prec * x
mvg_precmult(g::MvDG, x::Vector{AbstractFloat}) =  x ./ g.sigma2

# -----------------------------------------------------------------------------

"""
    gradloglik(g, w)

Gradient of the loglikelihood of a multivariate Gaussian `g` at point `x`.
"""
gradloglik(g::MvGaussian, x::Vector{AbstractFloat}) =
    mvg_precmu(g) - mvg_precmult(g, x)


"""
    loglik(g, x)

Loglikelihood of a multivariate Gaussian `g` at point `x`.
"""
loglik(g::MvGaussian, x::Vector{AbstractFloat}) =
    0.5dot(mvg_precmu(g) - mvg_precmult(g, x), x - mvg_mu(g))
