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
#> 1           random_error  0.003412853  0.003352479  0.001830899  0.005457652
#> 2          twofold_error  2.024316610  2.037772925  1.777694542  2.238514748
#> 3      boosting_children  2.725591888  2.748241417  1.854077602  3.484546466
#> 4        waning_children  0.120885888  0.115391366  0.059849537  0.201081699
#> 5        boosting_adults  2.488394977  2.483463123  2.021878424  2.981259657
#> 6          waning_adults  0.055891705  0.058434216  0.034719929  0.076565611
#> 7      inf_prob_children  0.542112263  0.542481295  0.367732314  0.745727722
#> 8        inf_prob_adults  0.354311038  0.356557796  0.272281582  0.429585140
#> 9  inf_prob_older_adults  0.254116503  0.248833975  0.166219692  0.383262753
#> 10              hai_coef -0.488309214 -0.492834293 -0.653790639 -0.323175360
# }
```
