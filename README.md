<!-- README.md is generated from README.Rmd. Please edit that file -->
tbn : A package to deal with Tweedie Bayesian networs
=====================================================

This package provides several useful functions to deal with Tweedie
Regression mode (TRMs) and Tweedie Bayesian Networks (TBNs).
$$ \\begin{equation} f(x\|\\Theta) = \\sum\_{k=1}^{K}{\\pi\_k f\_k(x\|\\mu\_k,\\Sigma\_k)}\\end{equation}$$

where *f*<sub>*k*</sub>(*x*\|*μ*<sub>*k*</sub>, *Σ*<sub>*k*</sub>) is
the density function of the multivariate normal distribution with mean
*μ*<sub>*k*</sub> and covariance matrix *Σ*<sub>*k*</sub>, and the
mixing proportions 0 &lt; *π*<sub>*k*</sub> &lt; 1 satisfy
$\\displaystyle\\sum\_{k=1}^{K}\\pi\_{k}=1$. In addition, each component
of this mixture is associated with a decomposable undirected graph
*G*<sub>*k*</sub> = (*V*, ℰ<sub>*k*</sub>), where *V* is the vertices
(nodes) set and ℰ<sub>*k*</sub> corresponds to the edges of the graph
*G*<sub>*k*</sub>. The set of all the mixture parameters is
*Θ* = {*π*<sub>1</sub>, ..., *π*<sub>*K*</sub>, *μ*<sub>1</sub>, ..., *μ*<sub>*K*</sub>, *Σ*<sub>1</sub>, ..., *Σ*<sub>*K*</sub>}

This package allows the estimation of the mixture parameters and data
classification. These tasks are achieved using an extended Expectation
Maximization algorithm called Graphical Expectation Maximization (GEM)
algorithm.

This package exports the following functions:

-   graphSigma
-   graphMatrixAssoc
-   computeTau
-   gemEstimator.

Required set-up for this package
--------------------------------

Currently, this package exists in a development version on GitHub. To
use the package, you need to install it directly from GitHub using the
`install_github` function from `devtools`.

You can use the following code to install the development version of
`countyweather`:

``` r
library(devtools)
install_github("km20/gemalogrithm")
library(gemalgorithm)
```

Applying Lauritzen’s formula : graphSigma
-----------------------------------------

This function applies the lauritzen’s formula to a covariance matrix to
take into account a decomposable graph structure. It uses the provided
covariance matrix and the provided graph to compute the covariance
matrix that perfectly fits the set of conditional independence
relationships encoded by the graph’s cliques and seperators.

### Example :

``` r
A <- matrix(0,5,5)
diag(A) <- 1:5
diag(A[-1,-5]) <- 1:4
diag(A[-5,-1]) <- 1:4
print(A)
cliques <- list(c(1,2),c(2,3), c(3,4),c(4,5))
separators <- list(c(2),c(3),c(4))
nS <- c(1,1,1)
Anew <- graphSigma(A, cliques, separators, nS)
print(Anew)
```

Association degree between an undirected graph and a covariance matrix: graphMatrixAssoc
----------------------------------------------------------------------------------------

This function computes the association degree between a covariabce
matrix and a graph. The computed metric relies only on the
correspondance between the zeros in the inverted covariance matrix and
the set of conditional independencies.

``` r
d1 <- graphMatrixAssoc(A,cliques)
d0 <- graphMatrixAssoc(Anew,cliques)
```

Since Anew prefectly matches the conditional independencies in the
graph, d0 is equal to 0. However, using the original matrix A, we get a
value of d1 equal to 1.95.

Posterior probability : computeTau
----------------------------------

The “computeTau” function calculates the posterior probability that each
observation belongs to each of the mixture components.
*τ*<sub>*i**j*</sub> is the posterior probability that the observation
*X*<sub>*i*</sub> belongs to the *j*<sup>*t**h*</sup> component of the
mixture and given by: $ \_{ij} = $

This function returns a matrix with n rows ( observations number) and K
columns (mixture components number).

Parameters estimation : gemEstimator
------------------------------------

The main function in this package is the “gemEstimator” which estimates
the Gaussian mixture parameters using the GEM algorithm. The mixture
components number is supposed to be known. This function uses an initial
parameters guess and a set of associated graphs to iteratively estimate
the parameters.

Starting from an intial parameters set *Θ*<sup>(0)</sup>, this function
repeats iteratively the 3 steps of the GEM algorithm :

-   Expectation step : Computes the conditional expectation of the
    complete-data log-likelihood given the observed data, using the
    current fit *Θ*<sup>(*l*)</sup> :
    *Q*(*Θ*\|\|*Θ*<sup>(*l*)</sup>) = *E*<sub>*Θ*<sup>(*l*)</sup></sub>(*L*(*X*<sub>1</sub>, ..., *X*<sub>*n*</sub>, *Z*<sub>1</sub>, ..., *Z*<sub>*n*</sub>, *Θ*)\|*X*<sub>1</sub>, ..., *X*<sub>*n*</sub>)

-   Maximization step: Consists in a global maximization of
    *Q*(*Θ*\|\|*Θ*<sup>(*l*)</sup>) with respect to *Θ* :
    *Θ*<sup>(*l* + 1)</sup> = arg max<sub>*Θ*</sub>*Q*(*Θ*\|\|*Θ*<sup>(*l*)</sup>)

-   G-Step : Applies the Lauriten’s formula to the estimated covariance
    matrices in order to take into account the known independencies.

The stopping rule depends on the “Nit” parameter used in the function
gemEstimator:

-   If Nit &gt; 0 : The algorithm stops after exactly Nit iterations.
-   If Nit &lt; 0 : The algorithm stops when :
    $$
    \\frac{\|\|\\Theta^{l+1} -\\Theta^{l}\|\|}{1+\|\|\\Theta^{l}\|\|} &lt; 10^{-4}
    $$
