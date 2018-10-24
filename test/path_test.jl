using PDSampler
using Random

xs = [[0.,0.] [1.,1.] [2.,1.5]]
ts = [0.0, 1.0, 2.0]

pa = Path(xs,ts)

@test samplepath(pa, 0.3) == [0.3, 0.3] &&
      samplepath(pa, 1.3) == [1.3, 1.15]

# Seg1: [0.0, 1.0]
# Seg2: [1.0, 2.0]

@test norm(samplepath(pa, [0.3, 1.3, 1.4]) -
            [0.3 1.3 1.4; 0.3 1.15 1.2]) <= 1e-12

@test_throws AssertionError samplepath(pa, 10.)

xs = [[0., 0.] [1., 1.] [2., 1.5] [3., 5.0] [3.4, -2.0]]
ts = [0.0, 1.0, 2.0, 2.5, 2.7]

pb = Path(xs, ts)

N  = 10000
ss = range(0.0, stop=ts[end], length=N)

samples = samplepath(pb, ss)
@test norm(sum(samples, dims=2) / N - pathmean(pb)) / norm(pathmean(pb)) <= 5e-3

#####
Random.seed!(15)
N  = 500
xs = zeros(2, N)
ts = zeros(N)
xs[:, 1] = randn(2)
for s ∈ 2:N
    xs[:, s] = randn(2)
    ts[s] = ts[s-1] + rand()
end
pb = Path(xs, ts)
N  = 10000
ss = range(0.0, stop=ts[end], length=N)
samples = samplepath(pb, ss)
@test norm(sum(samples, dims=2) / N - pathmean(pb)) / norm(pathmean(pb)) <= 5e-3
