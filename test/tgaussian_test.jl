using Random
using LinearAlgebra

# Everything should be the same as usual normals except
# when x is out of boundaries

p = 50
a = zeros(p)
b = ones(p)

mu  = rand(p)
P1  = randn(p, p)
P1 *= P1'
C1  = inv(P1)
C1 += C1'
C1 /= 2

mvg   = MvGaussianCanon(mu, P1)
mvg2  = MvGaussianNatural(P1*mu, -P1)
mvg_d = MvDiagonalGaussian(mu, 2ones(p))
mvgs  = MvGaussianStandard(mu, C1)

tmvg   = TMvGaussianCanon(mu, P1, a, b)
tmvg2  = TMvGaussianNatural(P1*mu, -P1, a, b)
tmvg_d = TMvDiagonalGaussian(mu, 2ones(p), a, b)
tmvgs  = TMvGaussianStandard(mu, C1, a, b)

x_valid = rand(p)

@test loglik(mvg, x_valid) == loglik(tmvg, x_valid)
@test loglik(mvg2, x_valid) == loglik(tmvg2, x_valid)
@test loglik(mvg_d, x_valid) == loglik(tmvg_d, x_valid)
@test loglik(mvgs, x_valid) == loglik(tmvgs, x_valid)

@test gradloglik(mvg, x_valid) == gradloglik(tmvg, x_valid)
@test gradloglik(mvg2, x_valid) == gradloglik(tmvg2, x_valid)
@test gradloglik(mvg_d, x_valid) == gradloglik(tmvg_d, x_valid)
@test gradloglik(mvgs, x_valid) == gradloglik(tmvgs, x_valid)

x_not = -rand(p)

@test loglik(tmvg, x_not) == -Inf
@test loglik(tmvg2, x_not) == -Inf
@test loglik(tmvg_d, x_not) == -Inf
@test loglik(tmvgs, x_not) == -Inf

infvec = Inf*ones(p)

@test gradloglik(tmvg, x_not) == infvec
@test gradloglik(tmvg2, x_not) == infvec
@test gradloglik(tmvg_d, x_not) == infvec
@test gradloglik(tmvgs, x_not) == infvec
