using PDSampler
using Base.Test

# the functions are somewhat trivial so the tests are also a bit trivial

Random.seed!(1240)

d = 10
normal = randn(d)
v1 = randn(d)
v2 = copy(v1)
I = eye(d)

@test norm(reflect_bps!(normal,v1)-reflect_bps!(normal,v2,I)) <= 1e-12

v3 = v1-2.0dot(normal,v1)*normal/norm(normal)^2

@test norm(reflect_bps!(normal,v1)-v3) <= 1e-12

mask     = rand(1:d, 3)
v4       = copy(v1)
v4[mask] = v4[mask]*(-1.0)

@test norm(reflect_zz!(mask,v1) - v4) <= 1e-12

v5          = copy(v1)
v5[mask[1]] = -v5[mask[1]]

@test norm(reflect_zz!(mask[1], v1) - v5) <= 1e-12

Random.seed!(12); v  = randn(d)
Random.seed!(12); v1 = refresh_global!(v1)

@test norm(v-v1) <= 1e-12

Random.seed!(31); v  = randn(d); v /= norm(v)
Random.seed!(31); v1 = refresh_restricted!(v1)

@test norm(v-v1) <= 1e-12

beta  = Beta(abs(randn()))
Random.seed!(53)
w  = randn(d)
w /= norm(w)
w *= tan( rand(beta) * 2 * pi )
v  = v1 + w
v /= norm(v)
Random.seed!(53)
v1 = refresh_partial!(beta)(v1)

@test norm(v-v1) <= 1e-12
