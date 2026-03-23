# Compute MCMC summary statistics

Compute MCMC summary statistics

## Usage

``` r
.mcmc_summary(mcmc_matrix)
```

## Arguments

- mcmc_matrix:

  Matrix of posterior samples (rows = iterations, cols = parameters).

## Value

A matrix with `ncol(mcmc_matrix)` rows and 4 columns: mean, 2.5%
quantile, 97.5% quantile, and acceptance rate.
