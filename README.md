
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mpower

<!-- badges: start -->

<!-- badges: end -->

The package mpower provides a framework for doing power analysis using
Monte Carlo simultions that

  - calculates statistical power for a wide range of models.

  - generates synthetic data with dependence structures of interests.

  - adjusts “effect size” using the signal-to-noise ratio that is
    applicable to continuous/discrete outcomes and non-linear
    conditional mean.

## Installation

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("phuchonguyen/mpower")
```

## Example

This is a basic example which shows you how to set up power analysis for
four correlated predictors, a continuous outcome, and using a linear
regression:

``` r
library(mpower)
# create four predictors with pairwise correlations around 0.6
xmod  <- MixtureModel(S=0.6, m=100, p=4, method="cvine",
                      var.name=c("x1", "x2", "x3", "x4"))
ymod  <- OutcomeModel(f="2*x1 + x2")
imod  <- InferenceModel(fun="glm")
# run power analysis with 1000 Monte Carlo iterations using two core
res <- sim_power(s=1000, n=100, sigma=1,
                xmod=xmod, ymod=ymod, imod=imod, cores=2)
```

We can now calculate the power from simulation samples:

``` r
summarize_power(res, crit = "ci")
#> For a 0.95 % CI
#>    x1    x2    x3    x4 
#> 1.000 1.000 0.052 0.046
```
