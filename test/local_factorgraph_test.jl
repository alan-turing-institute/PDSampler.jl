using PDMP, Base.Test

srand(1234)

p   = 2
P1  = randn(p)
P1 *= P1'
mu  = zeros(p)

mvg            = PDMP.MvGaussianCanon(mu,P1)
gll(x)         = PDMP.gradloglik(mvg, x)
nextevent(x,v) = PDMP.nextevent_bps(mvg, x, v)

chain = chaingraph([ Factor( nextevent, gll, i) for i in 1:3])


fgs = FactorGraphStruct( [[1,2,5],
                          [1,3,6],
                          [3,4,5,6],
                          [4,5]] )

@test assocfactors(fgs,3) == [2,3]
@test assocvariables(fgs,3) == [3,4,5,6]
@test sort(linkedfactors(fgs,3)) == [1,2,4]
@test sort(linkedfactors(fgs,4)) == [1,3]
@test chainstruct(5).flist == [[1,2],[2,3],[3,4],[4,5]]

xtest = randn(p)
vtest = randn(p)

srand(12); ne1 = chain.factors[3].nextevent(xtest,vtest)
srand(12); ne2 = nextevent(xtest,vtest)

@test   chain.structure.flist == [[1,2],[2,3],[3,4]] &&
        chain.factors[1].gll(xtest) == gll(xtest) &&
        ne1.tau == ne2.tau

chain = chainstruct(3)
@test assocvariables(chain, 1) == [1,2] && assocvariables(chain, 2) == [2,3]
