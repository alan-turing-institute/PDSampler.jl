# About PDMP samplers

This page aims at giving a very brief introduction to the concept of PDMP samplers and give insight of how it is implemented although we cover the implementation in more details in the technical documentation. This is not meant to be a rigorous presentation of the algorithm (for this, please see the references at the bottom of this page). Rather, we focus here on the "large building blocks" behind the algorithm.

The purpose of the algorithm is to be able to evaluate expected values with respect to an arbitrary target distribution which we assume admits a probability density function  $\pi$. For simplicity, we assume that $\pi:C\to \mathbb R^+$ with $C\subseteq \mathbb R^p$, convex. The objective is therefore to compute a weighted integral of the form:

\begin{equation}
    \mathbb E_{\pi}[\varphi(X)] = \int_{C} \varphi(x)\pi(x)\,\mathrm{d}x
\end{equation}





### Global Samplers


### Local Samplers

### References

* Alexandre Bouchard-Côté, Sebastian J. Vollmer and Arnaud Doucet, [*The Bouncy Particle Sampler: A Non-Reversible Rejection-Free Markov Chain Monte Carlo Method*](https://arxiv.org/abs/1510.02451), arXiv preprint, 2015.
* Joris Bierkens, Alexandre Bouchard-Côté, Arnaud Doucet, Andrew B. Duncan, Paul Fearnhead, Gareth Roberts and Sebastian J. Vollmer, [*Piecewise Deterministic Markov Processes for Scalable Monte Carlo on Restricted Domains*](https://arxiv.org/pdf/1701.04244.pdf), arXiv preprint, 2017.
* Joris Bierkens, Paul Fearnhead and Gareth Roberts, [*The Zig-Zag Process and Super-Efficient Sampling for Bayesian Analysis of Big Data*](https://arxiv.org/pdf/1607.03188.pdf), arXiv preprint, 2016.
