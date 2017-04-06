using PDMP, Base.Test

# ------------------------------------------------------------------------------
# Quick version of functions (for testing)
function q_linearupperbound(xstar,gllstar,b,x,v)
    a = max(dot(-gllstar,v),0.) + norm(x-xstar)*b
    tauprime = -a/b + sqrt((a/b)^2+2randexp()/b)
    (tauprime, a+b*tauprime)
end

# Benchmark dirty version (for testing)
function q_sampletaugaussian(mu,C,x,v)
    r      = randexp()
    vciv   = dot(v,C\v)
    vcixmu = dot(v,C\(x-mu))
    t      = 0
    if vcixmu < 0 ## can be optimised because we now apriori
        t= -vcixmu/vciv+sqrt(2*r/vciv)
    else
        t= -vcixmu/vciv+sqrt((vcixmu/vciv)^2+2*r/vciv)
    end
    t
end
# ------------------------------------------------------------------------------

srand(123)

p = 5

gllstar = randn(p)
xstar   = randn(p)
x       = randn(p)
v       = randn(p)
g       = randn(p)
b       = randn()

lb = LinearBound(gllstar,xstar,b)

## check that q_linearupperbound and nextevent_bps match.

srand(12); (a1,a2) = q_linearupperbound(gllstar,xstar,b,x,v)
acc = rand() < -dot(g,v)/a2
srand(12); bounce  = nextevent_bps(lb,x,v)

@test abs(a1-bounce.tau) <= 1e-12 && bounce.dobounce(g,v) == acc

## check that q_sampletaugaussian and nextevent_bps match.

srand(123)

p   = 50
C1  = randn(p,p)
C1 *= C1'
mu  = randn(p)

mvg = MvGaussianStandard(mu, C1)

x = randn(p)
v = randn(p)
g = randn(p)

srand(12); tauq   = q_sampletaugaussian(mu, C1, x, v)
srand(12); bounce = nextevent_bps(mvg, x, v)

@test abs(tauq - bounce.tau) <= 1e-12 && bounce.dobounce(g,v)
