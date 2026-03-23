# Validate inputs for sero_reconstruct

Validate inputs for sero_reconstruct

## Usage

``` r
.validate_inputs(
  inputdata,
  inputILI,
  n_iteration,
  burnin,
  thinning,
  group_by = NULL
)
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

- group_by:

  Optional formula; when non-NULL, age_group value checks are skipped.

## Value

The (possibly column-renamed) `inputdata`.
