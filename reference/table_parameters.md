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
#> 1           random_error  0.001897294  0.001835798  0.001206076  0.002943005
#> 2          twofold_error  1.967354047  1.971187224  1.797104254  2.168819614
#> 3      boosting_children  3.229963924  3.281385035  2.600976635  3.811797935
#> 4        waning_children  0.054187302  0.053413644  0.020274494  0.093470877
#> 5        boosting_adults  2.886245848  2.890986600  2.573971729  3.177674317
#> 6          waning_adults  0.144958376  0.148381843  0.052920826  0.220790151
#> 7      inf_prob_children  0.599507748  0.607457707  0.411553704  0.794285807
#> 8        inf_prob_adults  0.362219144  0.361593751  0.287874719  0.442789319
#> 9  inf_prob_older_adults  0.257929407  0.253914989  0.166606062  0.368575122
#> 10              hai_coef -0.486856121 -0.513136047 -0.689081070 -0.269118400
# }
```
