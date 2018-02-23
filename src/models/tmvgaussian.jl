export
    TMvGaussian,
    TMvGaussianStandard,
    TMvDiagonalGaussian,
    TMvGaussianCanon,
    TMvGaussianNatural

abstract type TMvGaussian end

"""
    TMvGaussianStandard(mu, cov, lbd, ubd)

Standard representation of a truncated multivariate Gaussian instantiated from
the mean and the covariance matrix as well as the boundaries.
Users should prefer the Canonical representation if they can use that.
"""
struct TMvGaussianStandard <: TMvGaussian
    mu::Vector{Real}   # Precision*mean
    cov::Matrix{Real}  # negative Precision
    lbd::Vector{Real}  # list of lower boundaries
    ubd::Vector{Real}  # list of upper boundaries
    p::Int  # Dimensions
    #
    prec::Matrix{Real}
    precmu::Vector{Real}
    function TMvGaussianStandard(mu, cov, lbd, ubd)
        @assert length(mu) == size(cov, 1) == size(cov, 2)
        @assert length(mu) == length(lbd) == length(ubd)
        @assert all(lbd .<= mu .<= ubd)
        L = cholfact(Symmetric(cov))[:L] # stored+considered as lower triangular
        Li = inv(L)
        prec = Li'*Li # going through CholFact is more stable + guarantees PD
        new(mu, cov, lbd, ubd, length(mu), prec, prec * mu)
    end
end

"""
    TMvDiagonalGaussian(mu, sigma, lbd, ubd)

Representation of a truncated multivariate diagonal gaussian (diagonal
covariance) as well as the boundaries.
"""
struct TMvDiagonalGaussian <: TMvGaussian
    mu::Vector{Real}     # mean
    sigma::Vector{Real}  # standard deviation
    sigma2::Vector{Real} # standard deviation squared
    lbd::Vector{Real}    # list of lower boundaries
    ubd::Vector{Real}    # list of upper boundaries
    p::Int # dimensions
    function TMvDiagonalGaussian(mu, sigma, lbd, ubd)
        @assert length(mu) == length(lbd) == length(ubd)
        @assert all(lbd .<= mu .<= ubd)
        new(mu, sigma, sigma.^2, lbd, ubd, length(mu))
    end
end

"""
    TMvGaussianCanon(mu, prec, lbd, ubd)

Canonical representation of a truncated multivariate Gaussian instantiated from
the mean and the precision matrix as well as the boundaries.
This is the preferred representation.
"""
struct TMvGaussianCanon <: TMvGaussian
    mu::Vector{Real}     # mean
    prec::Matrix{Real}   # Precision
    precmu::Vector{Real} # Precision*mean
    lbd::Vector{Real}    # list of lower boundaries
    ubd::Vector{Real}    # list of upper boundaries
    p::Int                # Dimensions
    function TMvGaussianCanon(mu, prec, lbd, ubd)
        @assert length(mu)==size(prec,1)==size(prec,2)
        @assert length(mu)==length(lbd)==length(ubd)
        @assert all(lbd .<= mu .<= ubd)
        new(mu, prec, prec * mu, lbd, ubd, length(mu))
    end
end

"""
    TMvGaussianNatural(precmu, negprec)

Natural parameter space representation of a truncated multivariate Gaussian.
"""
struct TMvGaussianNatural <: TMvGaussian
    precmu::Vector{Real}   # Precision*mean
    negprec::Matrix{Real}  # negative Precision
    lbd::Vector{Real}      # list of lower boundaries
    ubd::Vector{Real}      # list of upper boundaries
    p::Int                  # Dimensions
    function TMvGaussianNatural(precmu, negprec, lbd, ubd)
        @assert length(precmu) == size(negprec,1) == size(negprec,2)
        @assert length(precmu) == length(lbd) == length(ubd)
        @assert all(lbd .<= -negprec\precmu .<= ubd)
        new(precmu, negprec, lbd, ubd, length(precmu))
    end
end

# ------------------------------------------------------------------------------
### Helper functions for efficient computations
# the type aliases are not exposed

const TMvGS = TMvGaussianStandard
const TMvGC = TMvGaussianCanon
const TMvGN = TMvGaussianNatural
const TMvDG = TMvDiagonalGaussian

mvg_mu(g::TMvGC) = g.mu
mvg_mu(g::TMvGN) = -g.negprec \ g.precmu
mvg_mu(g::TMvGS) = g.mu
mvg_mu(g::TMvDG) = g.mu

mvg_precmu(g::TMvGC) = g.precmu
mvg_precmu(g::TMvGN) = g.precmu
mvg_precmu(g::TMvGS) = g.precmu
mvg_precmu(g::TMvDG) = g.mu ./ g.sigma2

mvg_precmult(g::TMvGC, x::Vector{<:Real}) =  g.prec * x
mvg_precmult(g::TMvGN, x::Vector{<:Real}) = -g.negprec * x
mvg_precmult(g::TMvGS, x::Vector{<:Real}) =  g.prec * x
mvg_precmult(g::TMvDG, x::Vector{<:Real}) =  x ./ g.sigma2

# ------------------------------------------------------------------------------

"""
    gradloglik(g, w)

Gradient of the loglikelihood of a truncated multivariate Gaussian `g` at
point `x`.
"""
function gradloglik(g::TMvGaussian, x::Vector{<:Real})
    if all(g.lbd .<= x .<= g.ubd)
        return mvg_precmu(g) - mvg_precmult(g,x)
    else
        return Inf * ones(g.p)
    end
end

"""
    loglik(g, x)

Loglikelihood of a truncated multivariate Gaussian `g` at point `x`.
"""
function loglik(g::TMvGaussian, x::Vector{<:Real})
    if all(g.lbd .<= x .<= g.ubd)
        return 0.5dot(mvg_precmu(g) - mvg_precmult(g, x), x - mvg_mu(g))
    else
        return -Inf
    end
end
