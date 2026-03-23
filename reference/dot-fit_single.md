# Run MCMC for a single fit

Core fitting function used by
[`sero_reconstruct`](https://timktsang.github.io/seroreconstruct/reference/sero_reconstruct.md)
for both standard 3-age-group fits and single-group sub-fits.

## Usage

``` r
.fit_single(inputdata, inputILI, n_iteration, burnin, thinning, n_groups = 3L)
```

## Arguments

- inputdata:

  Data frame of individual-level data.

- inputILI:

  Data frame or matrix of influenza activity.

- n_iteration:

  Number of MCMC iterations.

- burnin:

  Burn-in iterations.

- thinning:

  Thinning interval.

- n_groups:

  Number of effective age groups (3 for standard, 1 for single-group).

## Value

A `seroreconstruct_fit` object.
