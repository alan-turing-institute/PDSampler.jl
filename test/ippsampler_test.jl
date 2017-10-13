using PDSampler, Base.Test, QuadGK.quadgk

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

### test NextEvent capsule (trivial test for code cov)

tau = 0.1
dobounce = (g,v)->-dot(g,v)<1.0
flipindex = 4
ne = NextEvent(tau;dobounce=dobounce,flipindex=flipindex)
g  = randn(4)
v  = randn(4)
@test ne.tau==tau && ne.dobounce(g,v)==dobounce(g,v) &&
        ne.flipindex==flipindex
ne2 = NextEvent(tau)
@test ne2.tau==tau && ne2.dobounce(g,v) &&
        ne2.flipindex==-1
ne3 = NextEvent(tau;dobounce=dobounce)
@test ne3.tau==tau && ne3.dobounce(g,v)==dobounce(g,v) &&
        ne3.flipindex==-1

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

## additional tests

nedef = NextEvent(0.5)
@test   nedef.tau == 0.5 && nedef.dobounce(randn(),randn()) &&
        nedef.flipindex == -1

srand(555); a  = nextevent_zz(mvg, x, v)
srand(555)
u1   = mvg.prec*x-mvg.precmu
u2   = mvg.prec*v
taus = zeros(p)
for i in 1:p
    ai = u1[i] * v[i]
    bi = u2[i] * v[i]
    ci = ai./bi
    ei = max(0.0,ci)*ci
    taus[i] = -ci + sqrt(ei + 2.0randexp()/abs(bi))
end
tau, flipindex = findmin(taus)

@test a.tau == tau && a.flipindex == flipindex

## testing nextevent_bps for PMF

d = 3
r = 3
s = 0.1

pmfg = PMFGaussian(r, s, d)

n  = 100
xu = [randn(d) for i in 1:n]
xv = [randn(d) for i in 1:n]
wu = [randn(d) for i in 1:n]
wv = [randn(d) for i in 1:n]

srand(542)
tau = [nextevent_bps(pmfg, [xu[i];xv[i]], [wu[i];wv[i]]).tau for i in 1:n]

srand(542)
ri = [randexp() for i in 1:100]

# it must collide
p1(xu,xv,wu,wv) = PDSampler.Poly( [  dot(xu,xv)-pmfg.r,
                                (dot(xu,wv)+dot(xv,wu)),
                                dot(wu,wv) ])
p2(xu,xv,wu,wv) = PDSampler.Poly( [  (dot(xu,wv)+dot(xv,wu)),
                                2.0dot(wu,wv)])
E(xu,xv,wu,wv) = p1(xu,xv,wu,wv) * p2(xu,xv,wu,wv)

#χ(i) = t->max(0.0,E(xu[i],xv[i],wu[i],wv[i])(t))
χ(i) = t-> max(0.0, p1(xu[i],xv[i],wu[i],wv[i])(t) *
                    p2(xu[i],xv[i],wu[i],wv[i])(t) )

### testing roots
t0(i) = -0.5(dot(xu[i],wv[i])+dot(xv[i],wu[i]))/dot(wu[i],wv[i])
Δ(i)  = t0(i)^2-(dot(xu[i],xv[i])-pmfg.r)/dot(wu[i],wv[i])

@test maximum( (Δ(i)>0) ? maximum(
                    map( χ(i), t0(i)+sqrt(Δ(i))*[-1.0,0.0,1.0])) :
                                χ(i)(t0(i))
                for i in 1:n) <= 1e-12

tm(i) = t0(i) - sqrt(abs(Δ(i)))
tp(i) = t0(i) + sqrt(abs(Δ(i)))

@test maximum( quadgk(χ(i),0.0,tau[i])[1]-ri[i] for i in 1:n ) <= 1e-6
