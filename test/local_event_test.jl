using PDMP, Base.Test

srand(1234)

p     = 2
P1    = randn(p)
P1   *= P1'
mu    = zeros(p)
mvg   = MvGaussianCanon(mu,P1)
chain = chaingraph([ Factor( (x,v)->nextevent_bps(mvg, x, v),
                              x   ->gradloglik(mvg, x),
                              i) for i in 1:3])

srand(53)

# test_local_definegraph defines a chain.
nvars = chain.structure.nvars
x0    = randn(nvars)
v0    = randn(nvars)
v0   /= norm(v0)

all_evlist = AllEventList(Float, nvars)

for i in 1:nvars
    evi = Event(x0[i], v0[i], 0.0)
    pushevent!(all_evlist.evl[i], evi)
end

i = rand(1:nvars)
t = rand()

@test isapprox(samplelocalpath(all_evlist.evl[i], t), x0[i]+t*v0[i])

t1 = t+rand()
ev1 = Event(randn(), randn(), t1)
pushevent!(all_evlist.evl[i], ev1)
t2 = t1+rand()

@test isapprox( samplelocalpath(all_evlist.evl[i], t), x0[i]+t*v0[i])
@test isapprox( samplelocalpath(all_evlist.evl[i], t2), ev1.x+(t2-ev1.t)*ev1.v)

t3  = t2+rand()
ev2 = Event(randn(), randn(), t2)
pushevent!(all_evlist.evl[i], ev2)

ts  = sort(rand(100))
ts /= maximum(ts)
ts *= t3
ss  = samplelocalpath(all_evlist.evl[i],ts)
ses = [ ts[j]<t1?(x0[i]+ts[j]*v0[i]):
            (ts[j]<t2?(ev1.x+(ts[j]-ev1.t)*ev1.v):
                (ev2.x+(ts[j]-ev2.t)*ev2.v)) for j in 1:length(ts)]

# sample localpath (multiple between)
@test isapprox( norm(ss-ses)/length(ss) , 0.0 )
