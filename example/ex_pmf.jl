using PDMP

a = readdlm("data/pmf_testR.csv", ',', Int64)
R = a[:,1:3]

latentD = 4     # dimension of latent space
sigmaU  = 1.0
sigmaV  = 1.0

### there may be discrepancy with lines missing etc.
# -> use unique
nU = length(unique(R[:,1]))
nV = length(unique(R[:,2]))

# factors: create N factors for the users,
mvgU             = PDMP.MvDiagonalGaussian(zeros(latentD), sigmaU)
gllU(x)          = PDMP.gradloglik(mvgU, x)
nexteventU(x, v) = PDMP.nextevent_bps(mvgU, x, v)
factorU(k)       = Factor(nexteventU, gllU, k)

allfactors = [factorU(k) for k in 1:nU]

# factors: create M factors for the movies,
mvgV             = PDMP.MvDiagonalGaussian(zeros(latentD), sigmaV)
gllV(x)          = PDMP.gradloglik(mvgV, x)
nexteventV(x, v) = PDMP.nextevent_bps(mvgV, x, v)
factorV(k)       = Factor(nexteventV, gllV, k)

allV = [factorV(k) for k in 1:nV]

push!(allfactors, allV...)

# factors: create nz(R) factors for the ratings
maskU(x) = x[1:latentD]
maskV(x) = x[latentD+1:end]

d = Dict{Tuple{Int,Int},Int}()

# for k in 1:size(R,1)
#     i,j       = R[k,1:2]
#     d[(i,j)]  = k+nU+nV
#     fR(ui,vj) = PDMP.PMFGaussian(ui,vj,)
#
#
#     push!(allfactors, facRk)
# end

println("success")
