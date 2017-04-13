using PDMP, Base.Test

# testing empty structures
zfactor = emptyfactor()
@test   zfactor.nextevent(Inf)==Inf &&
        zfactor.gll(Inf)==Inf &&
        zfactor.index == 0

zfgs = emptyfactorgraphstruct()
@test isempty(zfgs.flist) && isempty(zfgs.vlist) &&
        zfgs.nfactors==0 && zfgs.nvars==0

zfg = emptyfactorgraph()
@test isempty(zfg.structure.flist) &&
        isempty(zfg.factors)

srand(1234)

p   = 2
P1  = randn(p)
P1 *= P1'
mu  = zeros(p)

mvg            = PDMP.MvGaussianCanon(mu,P1)
gll(x)         = PDMP.gradloglik(mvg, x)
nextevent(x,v) = PDMP.nextevent_bps(mvg, x, v)

listfactors = [ Factor( nextevent, gll, i) for i in 1:3]
chain       = chaingraph(listfactors)
compc       = FactorGraphStruct([[1,2],[2,3],[3,4]])
chain2      = FactorGraph(FactorGraphStruct([[1,2],[2,3],[3,4]]), listfactors)
chain3      = FactorGraph([[1,2],[2,3],[3,4]], listfactors)

i = rand(1:3)
x = randn(p)
v = randn(p)

@test   chain.structure.flist == compc.flist &&
        chain.structure.vlist == compc.vlist

@test   chain2.structure.flist == chain.structure.flist &&
        chain2.structure.vlist == chain.structure.vlist

@test   chain3.structure.flist == chain.structure.flist &&
        chain3.structure.vlist == chain.structure.vlist

@test chain.factors[i].gll(x) == chain2.factors[i].gll(x)

srand(12); a = chain2.factors[i].nextevent(x,v)
srand(12); b = chain.factors[i].nextevent(x,v)

@test a==b


fgs = FactorGraphStruct( [[1,2,5],
                          [1,3,6],
                          [3,4,5,6],
                          [4,5]] )

@test assocfactors(fgs,3) == [2,3]
@test assocvariables(fgs,3) == [3,4,5,6]
@test sort(linkedfactors(fgs,3)) == [1,2,4]
@test sort(linkedfactors(fgs,4)) == [1,3]
@test chainstruct(5).flist == [[1,2],[2,3],[3,4],[4,5]]

@test assocfactors(chain,  2) == [1,2]
@test linkedfactors(chain, 2) == [1,3]

xtest = randn(p)
vtest = randn(p)

srand(12); ne1 = chain.factors[3].nextevent(xtest,vtest)
srand(12); ne2 = nextevent(xtest,vtest)

@test   chain.structure.flist == [[1,2],[2,3],[3,4]] &&
        chain.factors[1].gll(xtest) == gll(xtest) &&
        ne1.tau == ne2.tau

chain = chainstruct(3)
@test assocvariables(chain, 1) == [1,2] && assocvariables(chain, 2) == [2,3]
