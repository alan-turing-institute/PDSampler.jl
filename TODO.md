This is a list complementing the opened issues for things that need to be improved or clarified (i.e. require a bit more than just a few lines of code).

## Algorithm

- Clarifying the stopping criterion (+ adding a clear discussion of the topic). In theory the stopping criterion should be the "duration" of the trajectory. In practice it may be easier to control the number of gradient evaluations or number of events etc. A discussion and consequent adaptation would be useful.

### Local BPS - specific

- There is/seems to be an equivalence between LBPS working well and some form of sparsity of the graphical model being present (i.e.: updates don't require actions on too many factors). This should be discussed in more details / clarified. In particular there may be cases for which BPS > LBPS and vice versa, a non-expert user may want to know this beforehand.

## Examples

- Reproduce the chain of Bivariate Gaussian example as in paper by Bouchard-Cote, Vollmer and Doucet and validate. It's just a chain with bivariate gaussians with precision `[1 a; a 1]` where `a` is between `0.1` and `0.9`.
- (other repo?) examples with comparison to STAN's NUTS sampler
- (currently ongoing) use LBPS, BPS for the probabilistic matrix factorisation problem. Compare with HMC, Gibbs. Gibbs should be completely out of the game and, in theory, HMC also. 
