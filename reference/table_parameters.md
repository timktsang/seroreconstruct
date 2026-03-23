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
#> MCMC complete in 24 seconds. Use summary() to view estimates.
table_parameters(fit)
#>                Parameter         Mean       Median        Lower        Upper
#> 1           random_error  0.002956828  0.002906974  0.001569037  0.004667191
#> 2          twofold_error  2.368352201  2.372144930  2.091025237  2.628865034
#> 3      boosting_children  2.079983276  1.971288290  1.291885093  3.008690179
#> 4        waning_children  0.130438323  0.125775530  0.090151560  0.213527271
#> 5        boosting_adults  2.363475991  2.343188056  2.105772584  2.662035013
#> 6          waning_adults  0.179781841  0.168442566  0.105646977  0.272219207
#> 7      inf_prob_children  0.533773311  0.545835397  0.378097094  0.685133753
#> 8        inf_prob_adults  0.346093016  0.347937399  0.271269592  0.413561070
#> 9  inf_prob_older_adults  0.203600966  0.199629898  0.128021401  0.304588933
#> 10              hai_coef -0.222374194 -0.231526875 -0.340605876 -0.085398501
# }
```
