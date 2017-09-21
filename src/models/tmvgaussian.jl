export
    TMvGaussian,
    TMvGaussianStandard,
    TMvDiagonalGaussian,
    TMvGaussianCanon,
    TMvGaussianNatural

@compat abstract type TMvGaussian end

"""
    TMvGaussianStandard(mu, cov, lbd, ubd)

Standard representation of a truncated multivariate Gaussian instantiated from
the mean and the covariance matrix as well as the boundaries.
Users should prefer the Canonical representation if they can use that.
"""
immutable TMvGaussianStandard <: TMvGaussian
    mu::Vector{Float}   # Precision*mean
    cov::Matrix{Float}  # negative Precision
    lbd::Vector{Float}  # list of lower boundaries
    ubd::Vector{Float}  # list of upper boundaries
    p::Int              # Dimensions
    #
    prec::Matrix{Float}
    precmu::Vector{Float}
    function TMvGaussianStandard(mu, cov, lbd, ubd)
        @assert length(mu)==size(cov,1)==size(cov,2)
        @assert length(mu)==length(lbd)==length(ubd)
        @assert all(lbd .<= mu .<= ubd)
        L    = cholfact(Symmetric(cov))[:L] # stored+considered as lower triangular
        Li   = inv(L)
        prec = Li'*Li # going through CholFact is more stable + guarantees PD
        new(mu,cov,lbd,ubd,length(mu),prec,prec*mu)
    end
end

"""
    TMvDiagonalGaussian(mu, sigma, lbd, ubd)

Representation of a truncated multivariate diagonal gaussian (diagonal
covariance) as well as the boundaries.
"""
immutable TMvDiagonalGaussian <: TMvGaussian
    mu::Vector{Float}     # mean
    sigma::Vector{Float}  # standard deviation
    sigma2::Vector{Float} # standard deviation squared
    lbd::Vector{Float}    # list of lower boundaries
    ubd::Vector{Float}    # list of upper boundaries
    p::Int                # dimensions
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
immutable TMvGaussianCanon <: TMvGaussian
    mu::Vector{Float}     # mean
    prec::Matrix{Float}   # Precision
    precmu::Vector{Float} # Precision*mean
    lbd::Vector{Float}    # list of lower boundaries
    ubd::Vector{Float}    # list of upper boundaries
    p::Int                # Dimensions
    function TMvGaussianCanon(mu, prec, lbd, ubd)
        @assert length(mu)==size(prec,1)==size(prec,2)
        @assert length(mu)==length(lbd)==length(ubd)
        @assert all(lbd .<= mu .<= ubd)
        new(mu, prec, prec*mu, lbd, ubd, length(mu))
    end
end

"""
    TMvGaussianNatural(precmu, negprec)

Natural parameter space representation of a truncated multivariate Gaussian.
"""
immutable TMvGaussianNatural <: TMvGaussian
    precmu::Vector{Float}   # Precision*mean
    negprec::Matrix{Float}  # negative Precision
    lbd::Vector{Float}      # list of lower boundaries
    ubd::Vector{Float}      # list of upper boundaries
    p::Int                  # Dimensions
    function TMvGaussianNatural(precmu, negprec, lbd, ubd)
        @assert length(precmu) == size(negprec,1) == size(negprec,2)
        @assert length(precmu) == length(lbd) == length(ubd)
        @assert all(lbd .<= -negprec\precmu .<= ubd)
        new(precmu,negprec,lbd,ubd,length(precmu))
    end
end

# ------------------------------------------------------------------------------
### Helper functions for efficient computations
# the type aliases are not exposed

const TMvGS = TMvGaussianStandard
const TMvGC = TMvGaussianCanon
const TMvGN = TMvGaussianNatural
const TMvDG = TMvDiagonalGaussian

mvg_mu(g::TMvGC)::Vector{Float} = g.mu
mvg_mu(g::TMvGN)::Vector{Float} = -g.negprec\g.precmu
mvg_mu(g::TMvGS)::Vector{Float} = g.mu
mvg_mu(g::TMvDG)::Vector{Float} = g.mu

mvg_precmu(g::TMvGC)::Vector{Float} = g.precmu
mvg_precmu(g::TMvGN)::Vector{Float} = g.precmu
mvg_precmu(g::TMvGS)::Vector{Float} = g.precmu
mvg_precmu(g::TMvDG)::Vector{Float} = g.mu./g.sigma2

mvg_precmult(g::TMvGC, x::Vector{Float})::Vector{Float} =  g.prec * x
mvg_precmult(g::TMvGN, x::Vector{Float})::Vector{Float} = -g.negprec * x
mvg_precmult(g::TMvGS, x::Vector{Float})::Vector{Float} =  g.prec * x
mvg_precmult(g::TMvDG, x::Vector{Float})::Vector{Float} =  x./g.sigma2

# ------------------------------------------------------------------------------

"""
    gradloglik(g, w)

Gradient of the loglikelihood of a truncated multivariate Gaussian `g` at
point `x`.
"""
function gradloglik(g::TMvGaussian, x::Vector{Float})::Vector{Float}
    if all(g.lbd .<= x .<= g.ubd)
        return mvg_precmu(g) - mvg_precmult(g,x)
    else
        return Inf*ones(g.p)
    end
end

"""
    loglik(g, x)

Loglikelihood of a truncated multivariate Gaussian `g` at point `x`.
"""
function loglik(g::TMvGaussian, x::Vector{Float})::Float
    if all(g.lbd .<= x .<= g.ubd)
        return 0.5dot(mvg_precmu(g)-mvg_precmult(g,x),x-mvg_mu(g))
    else
        return -Inf
    end
end
