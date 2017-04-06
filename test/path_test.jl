using PDMP, Base.Test

xs = [[0.,0.] [1.,1.] [2.,1.5]]
ts = [0.0, 1.0, 2.0]

p1 = PDMP.Path(xs,ts)

@test samplepath(p1,0.3) == [0.3,0.3 ] && samplepath(p1,1.3) == [1.3,1.15]

# Seg1: [0.0, 1.0]
# Seg2: [1.0, 2.0]

@test norm(samplepath(p1, [0.3, 1.3, 1.4])-[0.3 1.3 1.4; 0.3 1.15 1.2]) <= 1e-12

@test_throws AssertionError samplepath(p1,10.)
