
<!-- README.md is generated from README.Rmd. Please edit that file -->

# seroreconstruct

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

`seroreconstruct` is a Bayesian modeling framework to infer influenza
virus infection status, antibody dynamics, and individual infection
risks from serological data, by accounting for measurement error. This
could identifiy influenza infections by relaxing 4-fold rise rule, and
quantifies the contributions of age and pre-epidemic
hemagglutination-inhibiting (HAI) titers to infection risk.

While the package currently provides a set of fundamental functions, we
are actively working on expanding its capabilities with more advanced
tools for a comprehensive understanding of your results.

## Installation

1.  Install \[R\]\[r-project\]

2.  Install the development version of seroreconstruct from
    [GitHub](https://github.com/timktsang/seroreconstruct):

``` r
devtools::install_github("timktsang/seroreconstruct")
library(seroreconstruct)
```

## Example

This is a basic example of how to load some serological data and fitting
the model using the MCMC framework.

``` r
library(seroreconstruct)

## Load in a data set and flu activity data
data("inputdata")
data("flu_activity")

###### run the MCMC to estimate parameter of the model
###### in actual analysis, number of iteration is 200000, burnin is 100000, and thinning is 10
mcmc_result <- sero_reconstruct(inputdata,flu_activity, 2000,1000,1)

##### obtain the model estimate from fitted MCMC result
extract_mcmc_result <- summary(mcmc_result)
```

![The output of the MCMC results.](man/figures/sero_result.png)

## Citation

To cite package **seroreconstruct** in publications use:

Tsang TK, Perera RAPM, Fang VJ, Wong JY, Shiu EY, So HC, Ip DKM, Malik
Peiris JS, Leung GM, Cowling BJ, Cauchemez S. (2022). Reconstructing
antibody dynamics to estimate the risk of influenza virus infection. Nat
Commun. 2022 Mar 23;13(1):1557.
