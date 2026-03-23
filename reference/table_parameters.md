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
#> MCMC complete in 26 seconds. Use summary() to view estimates.
table_parameters(fit)
#>                Parameter         Mean       Median        Lower        Upper
#> 1           random_error  0.002854527  0.002861994  0.002101813  0.003769735
#> 2          twofold_error  2.257183635  2.260906134  2.065284603  2.467051963
#> 3      boosting_children  2.465443257  2.599924775  1.784319587  3.019179778
#> 4        waning_children  0.105325908  0.098678557  0.073117802  0.156954115
#> 5        boosting_adults  2.649053018  2.630014569  2.272943051  3.154906441
#> 6          waning_adults  0.225673563  0.216319565  0.136794318  0.357794969
#> 7      inf_prob_children  0.545233315  0.533351676  0.423580521  0.775363489
#> 8        inf_prob_adults  0.356582673  0.354190191  0.299828374  0.422715963
#> 9  inf_prob_older_adults  0.212112447  0.208986681  0.139227272  0.302212603
#> 10              hai_coef -0.308928834 -0.300653908 -0.471886798 -0.212798217
# }
```
