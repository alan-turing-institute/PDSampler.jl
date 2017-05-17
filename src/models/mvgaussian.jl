export
    MvGaussian,
    MvGaussianStandard,
    MvDiagonalGaussian,
    MvGaussianCanon,
    MvGaussianNatural

@compat abstract type MvGaussian end

"""
    MvGaussianStandard(mu, cov)

Standard representation of a multivariate Gaussian instantiated from the mean
and the covariance matrix. Users should prefer the Canonical representation if
they can use that.
"""
immutable MvGaussianStandard <: MvGaussian
    mu::Vector{Float}   # Precision*mean
    cov::Matrix{Float}  # negative Precision
    p::Int              # Dimensions
    #
    prec::Matrix{Float}
    precmu::Vector{Float}
    function MvGaussianStandard(mu, cov)
        L    = cholfact(cov)[:L] # stored+considered as lower triangular
        Li   = inv(L)
        prec = Li'*Li # going through CholFact is more stable + guarantees PD
        new(mu,cov,length(mu),prec,prec*mu)
    end
end

"""
    MvDiagonalGaussian(mu,sigma)

Representation of a multivariate diagonal gaussian (diagonal covariance).
"""
immutable MvDiagonalGaussian <: MvGaussian
    mu::Vector{Float} # mean
    sigma::Float      # standard deviation
    sigma2::Float     # standard deviation squared
    p::Int            # dimensions
    MvDiagonalGaussian(mu, sigma) = new(mu, sigma, sigma^2, length(mu))
end

"""
    MvGaussianCanon

Canonical representation of a multivariate Gaussian instantiated from the mean
and the precision matrix. This is the preferred representation.
"""
immutable MvGaussianCanon <: MvGaussian
    mu::Vector{Float}     # mean
    precmu::Vector{Float} # Precision*mean
    prec::Matrix{Float}   # Precision
    p::Int                # Dimensions
    function MvGaussianCanon(mu, prec)
        @assert length(mu)==size(prec,1)==size(prec,2)
        new(mu,prec*mu,prec,length(mu))
    end
end

"""
    MvGaussianNatural

Natural parameter space representation of a multivariate Gaussian.
"""
immutable MvGaussianNatural <: MvGaussian
    precmu::Vector{Float}   # Precision*mean
    negprec::Matrix{Float}  # negative Precision
    p::Int                  # Dimensions
    MvGaussianNatural(precmu,negprec) = new(precmu,negprec,length(precmu))
end

# ------------------------------------------------------------------------------
### Helper functions for efficient computations
# the type aliases are not exposed

const MvGS = MvGaussianStandard
const MvGC = MvGaussianCanon
const MvGN = MvGaussianNatural
const MvDG = MvDiagonalGaussian

mvg_mu(g::MvGC)::Vector{Float} = g.mu
mvg_mu(g::MvGN)::Vector{Float} = -g.negprec\g.precmu
mvg_mu(g::MvGS)::Vector{Float} = g.mu
mvg_mu(g::MvDG)::Vector{Float} = g.mu

mvg_precmu(g::MvGC)::Vector{Float} = g.precmu
mvg_precmu(g::MvGN)::Vector{Float} = g.precmu
mvg_precmu(g::MvGS)::Vector{Float} = g.precmu
mvg_precmu(g::MvDG)::Vector{Float} = g.mu/g.sigma2

mvg_precmult(g::MvGC, x::Vector{Float})::Vector{Float} =  g.prec * x
mvg_precmult(g::MvGN, x::Vector{Float})::Vector{Float} = -g.negprec * x
mvg_precmult(g::MvGS, x::Vector{Float})::Vector{Float} =  g.prec * x
mvg_precmult(g::MvDG, x::Vector{Float})::Vector{Float} =  x/g.sigma2

# ------------------------------------------------------------------------------

"""
    gradloglik(g, w)

Gradient of the loglikelihood of a multivariate Gaussian `g` at point `x`.
"""
function gradloglik(g::MvGaussian, x::Vector{Float})::Vector{Float}
    mvg_precmu(g) - mvg_precmult(g,x)
end

"""
    loglik(g, x)

Loglikelihood of a multivariate Gaussian `g` at point `x`.
"""
function loglik(g::MvGaussian, x::Vector{Float})::Float
    0.5dot(mvg_precmu(g)-mvg_precmult(g,x),x-mvg_mu(g))
end
