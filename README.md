<!-- README.md is generated from README.Rmd. Please edit that file -->
tbn : A package to deal with Tweedie Bayesian networks
======================================================

General description
-------------------

A collection of algorithms and functions to aid statistical modeling.
This package provides several useful functions to deal with both the
Tweedie Regression models (TRMs) and the Tweedie Bayesian Networks
(TBNs). (TBNs) constitutes a new family of continuous Bayesian networks
whose conditional distributions belong to the Tweedie class.

This package provides a full procedure for the learning of the TBN
structure and parameters. It also introduces a sensitivity measure
relying on the Kulback-Liebler divergence. The proposed algorithms
partially relies on functions from statmod and tweedie packages.

This package exports the following functions:

-   learn.tbn
-   learn.trm
-   scores.tbn
-   scores.trm
-   sensi.tbn
-   sensi.trm
-   generateSample

Required set-up for this package
--------------------------------

Currently, this package exists in a development version on GitHub. To
use the package, you need to install it directly from GitHub using the
`install_github` function from `devtools`.

You can use the following code to install the development version of
`countyweather`:

``` r
library(devtools)
install_github("km20/tbn")
library(tbn)
```
