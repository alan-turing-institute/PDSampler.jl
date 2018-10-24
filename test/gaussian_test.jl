# --------------------------------------------------------------------------
# Quick version of functions (for testing)
q_gradloglik(mu,cov,x) = cov\(mu-x)
# --------------------------------------------------------------------------

Random.seed!(123)

# Parameters for the gaussians
p   = 50         # dimensions
P1  = randn(p,p) # precision
P1 *= P1'        # make the precision positive def
mu  = randn(p)   # mean
C1  = inv(P1)    # covariance
C1 += C1'
C1 /= 2          # stabilise the covariance

# Define gaussians using the different notations
mvg       = MvGaussianCanon(mu, P1)
mvg2      = MvGaussianNatural(P1*mu,-P1)
mvg_distr = Distributions.MvNormal(mu, C1)
mvg_da    = MvGaussianCanon(mu, 0.25 * diagm(0=>ones(length(mu))))
mvg_db    = MvDiagonalGaussian(mu, 2.0*ones(p))

# we don't care about constants, but in Distributions
# they do so we need to compensate for that
check_loglik(x) = Distributions.logpdf(mvg_distr, x) +
                    0.5 * (p * log(2pi) + log(det(C1)))

xvals = [randn(p) for i in 1:25]

# difference is in logspace, and we're using slightly different matrices and
# the computation of the constant is crude to say the least -> higher tol
@test maximum(abs.(collect(loglik(mvg, x) - check_loglik(x) for x ∈ xvals))) <= 1e-6

# canon to natural
@test maximum(abs.(collect(loglik(mvg, x) - loglik(mvg2, x) for x ∈ xvals))) <= 1e-10

# standard to diagonal
@test maximum(abs.(collect(loglik(mvg_da, x) - loglik(mvg_db, x) for x ∈ xvals))) <= 1e-10

# gradloglik, expect a difference since C1 was stabilised
@test maximum(norm(gradloglik(mvg, x) - q_gradloglik(mu, C1, x))
                    for x in xvals ) <= 1e-6

# gradloglik standard to diagonal
@test maximum(norm(gradloglik(mvg_da, x) - gradloglik(mvg_db, x))
                    for x in xvals ) <= 1e-10
