This is a list complementing the opened issues for things that need to be improved or clarified (i.e. require a bit more than just a few lines of code).

## [C] Code base

- [C1] get rid of dependency on Klara.jl
- [C2] get rid of dependency on ApproxFuns.jl (or make it conditional in some way since very seldom use and slows down loading of PDMP hugely)
- [C3] get help from someone to use Compat properly so that the syntax is forward looking whilst not breaking on Julia 0.5 ideally (or breaking on 0.5 and working with 0.7)

## [A] Algorithm / Theory

- [A1] Clarifying the stopping criterion (+ adding a clear discussion of the topic). In theory the stopping criterion should be the "duration" of the trajectory. In practice it may be easier to control the number of gradient evaluations or number of events etc. A discussion and consequent adaptation would be useful.
- [A2] Clarify useful comparison rules with respect to other algorithms (e.g.: ESS/CPUsec).

### Local BPS - specific

- [A3] There is/seems to be an equivalence between LBPS working well and some form of sparsity of the graphical model being present (i.e.: updates don't require actions on too many factors). This should be discussed in more details / clarified. In particular there may be cases for which BPS > LBPS and vice versa, a non-expert user may want to know this beforehand.

## [D] Documentation

- [D1] Finish the discussion on sampling methods for general IPP (detail cases + examples)
- [D2] Discuss in details stopping criterions for BPS/LBPS (see also [A1])

## [E] Examples

- [E1] Reproduce the chain of Bivariate Gaussian example as in paper by Bouchard-Cote, Vollmer and Doucet and validate. It's just a chain with bivariate gaussians with precision `[1 a; a 1]` where `a` is between `0.1` and `0.9`.
- [E2] (other repo?) examples with comparison to STAN's NUTS sampler
- [E3] (currently ongoing) use LBPS, BPS for the probabilistic matrix factorisation problem. Compare with HMC, Gibbs. Gibbs should be completely out of the game and, in theory, HMC also.
