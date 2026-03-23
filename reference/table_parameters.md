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
#> 1           random_error  0.002365497  0.002270652  0.001648674  0.003554406
#> 2          twofold_error  2.065970006  2.064283722  1.864899582  2.252221437
#> 3      boosting_children  3.079938967  3.123297680  2.391657640  3.700671952
#> 4        waning_children  0.059152525  0.059571029  0.027047137  0.078111185
#> 5        boosting_adults  2.909670550  2.911762152  2.534913126  3.283571527
#> 6          waning_adults  0.209037074  0.213473568  0.117917849  0.325360681
#> 7      inf_prob_children  0.517863595  0.513696587  0.406800951  0.656930217
#> 8        inf_prob_adults  0.337942524  0.335058158  0.278445401  0.407005039
#> 9  inf_prob_older_adults  0.230724186  0.227079700  0.151263270  0.322479125
#> 10              hai_coef -0.356673833 -0.358703459 -0.469014812 -0.222111748
# }
```
