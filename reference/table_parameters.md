# Summary table of model parameters with credible intervals

Extracts posterior summaries (mean, median, credible intervals) for all
active model parameters.

## Usage

``` r
table_parameters(fit, probs = c(0.025, 0.975))
```

## Arguments

- fit:

  A `seroreconstruct_fit`, `seroreconstruct_joint`, or
  `seroreconstruct_multi` object.

- probs:

  Numeric vector of length 2 giving the lower and upper quantile
  probabilities for the credible interval. Default `c(0.025, 0.975)` for
  a 95% interval.

## Value

A data frame with columns: Parameter, Mean, Median, Lower, Upper.

## Examples

``` r
# \donttest{
fit <- sero_reconstruct(inputdata, flu_activity,
                        n_iteration = 2000, burnin = 1000, thinning = 1)
#> 1000
#> MCMC complete in 28 seconds. Use summary() to view estimates.
table_parameters(fit)
#>                Parameter         Mean       Median        Lower        Upper
#> 1           random_error  0.002223158  0.002144885  0.001445509  0.003227777
#> 2          twofold_error  2.113138595  2.109557613  1.934144709  2.342849078
#> 3      boosting_children  3.115977344  3.110639859  2.782389890  3.483439575
#> 4        waning_children  0.098018478  0.096173453  0.059914999  0.142621784
#> 5        boosting_adults  2.533081277  2.596249741  1.952154268  2.895433940
#> 6          waning_adults  0.180013040  0.185702942  0.119479332  0.231707432
#> 7      inf_prob_children  0.490786481  0.485409698  0.371360053  0.630993105
#> 8        inf_prob_adults  0.358487142  0.358872456  0.289750695  0.426806000
#> 9  inf_prob_older_adults  0.259958111  0.248000829  0.154912513  0.387684986
#> 10              hai_coef -0.348902295 -0.344081872 -0.438665487 -0.246203847
# }
```
