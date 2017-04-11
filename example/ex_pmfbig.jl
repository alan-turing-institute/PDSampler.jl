using PDMP, JLD

start = time()

a = readdlm("data/ratings_filtered.csv", ',', Int64)
R = a[:,1:3]

println("($(time()-start)s) -- read the data")

latentD = 50     # dimension of latent space
sigmaU  = 1.0
sigmaV  = 1.0
sigmaR  = 1.0

### there may be discrepancy with lines missing etc.
# -> use unique
nU = maximum(R[:,1])
nV = maximum(R[:,2])

# factors: create N factors for the users,
mvgU             = MvDiagonalGaussian(zeros(latentD), sigmaU)
gllU(x)          = gradloglik(mvgU, x)
nexteventU(x, v) = nextevent_bps(mvgU, x, v)
factorU(k)       = Factor(nexteventU, gllU, k)

allfactors = [factorU(k) for k in 1:nU]

# factors: create M factors for the movies,
mvgV             = MvDiagonalGaussian(zeros(latentD), sigmaV)
gllV(x)          = gradloglik(mvgV, x)
nexteventV(x, v) = nextevent_bps(mvgV, x, v)
factorV(k)       = Factor(nexteventV, gllV, nU+k)

allV = [factorV(k) for k in 1:nV]

push!(allfactors, allV...)

# factors: create nz(R) factors for the ratings
maskU(x) = x[1:latentD]
maskV(x) = x[latentD+1:end]

structure = [[k] for k in 1:(nU+nV)] # each factor connected to its own var

# -----------------------------
# factors  variables
#  fU1       U1 = [1]
#  ...       ...
#  fUnU      UnU = [nU]
#  fV1       V1  = [nU+1]
#  ...       ...
#  fVnV      VnV = [1+nV]
# (for appropriate ij)
#  fRij      Ui,Vj = [i, nU+j]
# -----------------------------

for k in 1:size(R,1)
    i,j, rij = R[k,:]
    # the likelihood
    gij = PMFGaussian(rij, sigmaR, latentD)
    # the factor
    fij = Factor( (x,w)->nextevent_bps(gij,x,w),
                   x->gradloglik(gij, x),
                   nU+nV+k )
    push!(allfactors, fij)
    push!(structure, [i, j+nU])
end

fg = FactorGraph(structure, allfactors)

lambdaref  = .01
maxnevents = 1000
T          = Inf
nvars      = nU+nV

x0 = [randn(latentD) for i in 1:nvars]
v0 = [randn(latentD) for i in 1:nvars]
v0 = map(v->v/norm(v), v0)

lsim = LocalSimulation(fg, x0, v0, T, maxnevents, lambdaref)

println("($(time()-start)s) -- created the graph + sim")

(all_evlist, details) = simulate(lsim)

println("($(time()-start)s) -- finished the simulation")

save("res/dex_lbps_pmfbig.jld", "evlist", all_evlist.evl, "details", details)

println("($(time()-start)s) -- saved the results")
