using PDMP, Base.Test

xs = [[0.,0.] [1.,1.] [2.,1.5]]
ts = [0.0, 1.0, 2.0]

pa = Path(xs,ts)

@test samplepath(pa,0.3) == [0.3,0.3 ] && samplepath(pa,1.3) == [1.3,1.15]

# Seg1: [0.0, 1.0]
# Seg2: [1.0, 2.0]

@test norm(samplepath(pa, [0.3, 1.3, 1.4])-[0.3 1.3 1.4; 0.3 1.15 1.2]) <= 1e-12

@test_throws AssertionError samplepath(pa,10.)

xs = [[0.,0.] [1.,1.] [2.,1.5] [3.,5.0] [3.4,-2.0]]
ts = [0.0, 1.0, 2.0, 2.5, 2.7]

pb = Path(xs,ts)

N  = 1000
T  = 2.6
ss = linspace(0.0,T,N)

samples = samplepath(pb,ss)
@test norm(sum(samples,2)/N - pathmean(pb,T)) <= 1e-3
