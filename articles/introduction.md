# Getting Started with seroreconstruct

## Overview

`seroreconstruct` is an R package for Bayesian inference of influenza
virus infection status from serial antibody measurements. It relaxes the
traditional 4-fold rise rule by jointly modeling antibody boosting after
infection, waning over time, and measurement error — estimating
individual-level infection probabilities from hemagglutination
inhibition (HAI) titer data.

The statistical framework is described in [Tsang et al. (2022) *Nature
Communications* 13:1557](https://doi.org/10.1038/s41467-022-29310-8).

## Input data

The package requires two data frames:

1.  **`inputdata`** — one row per individual with columns:
    - `age_group`: integer (0 = children, 1 = adults, 2 = older adults)
    - `start_time`, `end_time`: follow-up period (integer days)
    - `time1`, `time2`, `time3`: serum collection dates (integer days)
    - `HAI_titer_1`, `HAI_titer_2`, `HAI_titer3`: HAI titers at each
      collection
2.  **`inputILI`** — daily influenza-like illness activity, with rows
    corresponding to consecutive days. The row indices must cover the
    range of `start_time` to `end_time` in `inputdata`.

``` r
library(seroreconstruct)

data("inputdata")
data("flu_activity")

head(inputdata)
#>   age_group start_time end_time time1 time2 time3 HAI_titer_1 HAI_titer_2
#> 1         0        837     1192   837  1032  1192           0           0
#> 2         2        837     1192   837  1032  1192           0           0
#> 3         2        837     1192   837  1032  1192           0           0
#> 4         0        837     1192   837  1032  1192           4           4
#> 5         0        845     1191   845  1034  1191           0           0
#> 6         2        845     1191   845  1034  1191           0           0
#>   HAI_titer3
#> 1          0
#> 2          0
#> 3          0
#> 4          4
#> 5          0
#> 6          0
dim(inputdata)
#> [1] 1753    9
```

## Fitting the model

The main function
[`sero_reconstruct()`](https://timktsang.github.io/seroreconstruct/reference/sero_reconstruct.md)
runs a Bayesian MCMC to estimate infection probabilities and antibody
dynamics parameters. For real analyses, use at least 100,000 iterations
with appropriate burn-in and thinning.

``` r
fit <- sero_reconstruct(inputdata, flu_activity,
                        n_iteration = 200000, burnin = 100000, thinning = 10)
```

For this vignette, we use a short run for illustration:

``` r
fit <- sero_reconstruct(inputdata, flu_activity,
                        n_iteration = 5000, burnin = 2000, thinning = 1)
#> 1000
#> 2000
#> 3000
#> 4000
#> MCMC complete in 66 seconds. Use summary() to view estimates.
fit
#> seroreconstruct fit
#>   Individuals: 1753 
#>   Age groups: 3 
#>   Posterior samples: 3000 
#>   Runtime: 66 seconds
#> 
#> Use summary() to extract model estimates.
```

## Viewing results

The [`summary()`](https://rdrr.io/r/base/summary.html) method extracts
key estimates with 95% credible intervals:

``` r
summary(fit)
#>                                                                                                        Variable
#>                                                                                                Random error (%)
#>                                                                                              Two-fold error (%)
#>                                                           Fold-increase after infection for children (Boosting)
#>                                                                Fold-decrease after 1 year for children (Waning)
#>                                                                Fold-increase after 1 year for adults (Boosting)
#>                                                                  Fold-decrease after 1 year for adults (Waning)
#>                                                                              Infection probability for children
#>                                                                                Infection probability for adults
#>                                                                          Infection probability for older adults
#>                                             Infection probability for children with pre-epidemic HAI titer < 10
#>                                               Infection probability for adults with pre-epidemic HAI titer < 10
#>                                         Infection probability for older adults with pre-epidemic HAI titer < 10
#>                                                                        Relative risk for children (Ref: Adults)
#>                                                                    Relative risk for older adults (Ref: Adults)
#>      Relative risk for children with pre-epidemic HAI titer < 10 (Ref: Adults with pre-epidemic HAI titer < 10)
#>  Relative risk for older adults with pre-epidemic HAI titer < 10 (Ref: Adults with pre-epidemic HAI titer < 10)
#>  Point estimate Lower bound Upper bound
#>            1.57        1.04        2.60
#>            3.38        4.29        2.68
#>           10.89        7.45       15.39
#>            1.05        1.02        1.14
#>            7.12        4.83        9.61
#>            1.07        1.02        1.15
#>            0.19        0.16        0.23
#>            0.20        0.16        0.24
#>            0.13        0.08        0.19
#>            0.47        0.38        0.55
#>            0.30        0.25        0.35
#>            0.20        0.14        0.28
#>            0.98        0.76        1.23
#>            1.57        1.25        1.94
#>            0.67        0.43        0.96
#>            0.68        0.44        0.96
```

The summary table includes:

- **Random error (%)** and **Two-fold error (%)**: measurement error
  parameters
- **Boosting**: fold-increase in antibody titer after infection
- **Waning**: fold-decrease in antibody titer per year post-infection
- **Infection probability**: overall and stratified by baseline HAI
  titer

## Visualization

### MCMC diagnostics

Check convergence with trace plots and posterior density plots:

``` r
plot_diagnostics(fit, params = c("random_error", "twofold_error"))
```

![](introduction_files/figure-html/diagnostics-1.png)

### Antibody trajectories

Visualize posterior trajectories for individual participants. Red lines
show trajectories where infection occurred; blue lines show trajectories
without infection.

``` r
plot_trajectory(fit, id = 1)
```

![](introduction_files/figure-html/trajectory-1.png)

### Boosting and waning

``` r
par(mfrow = c(1, 2))
plot_boosting(fit)
plot_waning(fit)
```

![](introduction_files/figure-html/boosting-waning-1.png)

### Infection probabilities

Forest plot of posterior infection probabilities with 95% credible
intervals:

``` r
plot_infection_prob(fit)
```

![](introduction_files/figure-html/infection-prob-1.png)

## Tables

Extract parameter estimates and individual-level results as data frames:

``` r
# Model parameter estimates
table_parameters(fit)
#>                Parameter        Mean       Median        Lower        Upper
#> 1           random_error  0.00156971  0.001505829  0.001044175  0.002598653
#> 2          twofold_error  1.99958683  1.999856632  1.761938807  2.231197323
#> 3      boosting_children  3.44460186  3.453197949  2.896825461  3.944256450
#> 4        waning_children  0.07655252  0.064369825  0.032148707  0.194253834
#> 5        boosting_adults  2.83095125  2.859607355  2.271848792  3.265140615
#> 6          waning_adults  0.09836363  0.099350093  0.030406325  0.201842372
#> 7      inf_prob_children  0.68642876  0.683354967  0.513957120  0.868202810
#> 8        inf_prob_adults  0.38572199  0.383994061  0.312757064  0.472097303
#> 9  inf_prob_older_adults  0.24680940  0.245946594  0.158182661  0.350112895
#> 10              hai_coef -0.54061944 -0.543337155 -0.745716802 -0.347277526

# Per-individual infection probabilities
head(table_infections(fit))
#>   Individual Infection_prob Infection_time_mean Baseline_titer_mean
#> 1          1         0.0000                  NA                0.52
#> 2          2         0.0037               919.5                0.52
#> 3          3         0.0063               949.1                0.48
#> 4          4         0.0000                  NA                4.56
#> 5          5         0.0257               951.8                0.48
#> 6          6         0.0000                  NA                0.48
```

## Subgroup analysis

Compare infection rates across groups by fitting independent MCMCs:

``` r
fit_by_age <- sero_reconstruct(inputdata, flu_activity,
                               n_iteration = 200000, burnin = 100000,
                               thinning = 10, group_by = ~age_group)

summary(fit_by_age)
```

## Joint model with shared parameters

When comparing groups, measurement error is a lab assay property — it is
shared across all groups. You can additionally share boosting/waning if
the groups are expected to have similar antibody dynamics:

``` r
# Share all parameters except infection probability
fit_joint <- sero_reconstruct(inputdata, flu_activity,
                              n_iteration = 200000, burnin = 100000,
                              thinning = 10,
                              group_by = ~age_group,
                              shared = c("error", "boosting_waning"))

summary(fit_joint)
```

## Simulation

Generate synthetic data for power analysis or validation:

``` r
data("para1")
data("para2")

sim <- simulate_data(inputdata, flu_activity, para1, para2)
names(sim)
#> [1] "age_group"   "start_time"  "end_time"    "time1"       "time2"      
#> [6] "time3"       "HAI_titer_1" "HAI_titer_2" "HAI_titer3"
```

The simulated data can be passed directly to
[`sero_reconstruct()`](https://timktsang.github.io/seroreconstruct/reference/sero_reconstruct.md)
for simulation-recovery studies.

## Citation

To cite seroreconstruct in publications:

> Tsang TK, Perera RAPM, Fang VJ, Wong JY, Shiu EY, So HC, Ip DKM, Malik
> Peiris JS, Leung GM, Cowling BJ, Cauchemez S. (2022). Reconstructing
> antibody dynamics to estimate the risk of influenza virus infection.
> *Nature Communications* 13:1557. doi:
> [10.1038/s41467-022-29310-8](https://doi.org/10.1038/s41467-022-29310-8)
