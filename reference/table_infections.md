# Per-individual infection estimates

Summarizes posterior infection status, timing, and baseline titer for
each individual in the dataset.

## Usage

``` r
table_infections(fit)
```

## Arguments

- fit:

  A `seroreconstruct_fit` or `seroreconstruct_joint` object.

## Value

A data frame with one row per individual and columns: `Individual` (row
index), `Infection_prob` (posterior mean probability of infection),
`Infection_time_mean` (mean infection time among infected samples),
`Baseline_titer_mean` (mean imputed baseline HAI titer).

## Examples

``` r
# \donttest{
fit <- sero_reconstruct(inputdata, flu_activity,
                        n_iteration = 2000, burnin = 1000, thinning = 1)
#> 1000
#> MCMC complete in 26 seconds. Use summary() to view estimates.
head(table_infections(fit))
#>   Individual Infection_prob Infection_time_mean Baseline_titer_mean
#> 1          1           0.14               948.1                0.44
#> 2          2           0.00                  NA                0.51
#> 3          3           0.00                  NA                0.47
#> 4          4           0.00                  NA                4.54
#> 5          5           0.00                  NA                0.54
#> 6          6           0.00                  NA                0.49
# }
```
