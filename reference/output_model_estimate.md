# Extract the model estimates from the fitted MCMC

`output_model_estimate` is deprecated; use
[`summary()`](https://rdrr.io/r/base/summary.html) instead.

## Usage

``` r
output_model_estimate(fitted_MCMC, period)
```

## Arguments

- fitted_MCMC:

  A `seroreconstruct_fit` object, or a list returned by an older version
  of
  [`sero_reconstruct()`](https://timktsang.github.io/seroreconstruct/reference/sero_reconstruct.md).

- period:

  A vector indicating the start and the end of a season to compute the
  infection probabilities. If empty, the start and end of the season are
  inferred from the data.

## Value

A data frame of model estimates (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
a1 <- sero_reconstruct(inputdata, flu_activity,
                        n_iteration = 2000, burnin = 1000, thinning = 1)
fitted_result <- output_model_estimate(a1)  # deprecated, use summary(a1)
} # }
```
