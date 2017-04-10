## algorithm

- controlling the algorithm (stopping criterion matters, use Time)
- sparsity implied for efficiency of LBPS

## examples

- chain of bivariate gaussian with prec `[1 a; a 1]` where `a` is between `0.1` and `0.9` (paper)
- compare with STAN HMC

## PMF ratings (movielens etc)

Change of delimiters

```bash
a=ratings.dat; sed s/::/,/g $a > _tmp; mv _tmp $a
```

```julia
a = readdlm("ratings.dat",',',Int64)
```
