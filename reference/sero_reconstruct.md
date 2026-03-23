# Run the MCMC for the Bayesian model

The main function to run the MCMC for the Bayesian model, to obtain
individual dynamics, model parameters such as infection probability,
boosting, waning, and measurement error.

## Usage

``` r
sero_reconstruct(
  inputdata,
  inputILI,
  n_iteration = 2000,
  burnin = 1000,
  thinning = 1,
  group_by = NULL,
  shared = NULL,
  subject_ids = NULL
)
```

## Arguments

- inputdata:

  The data for running MCMC, in dataframe format. It should be in the
  same format as the data in the package. It includes: 1) age_group (0:
  children, 1: adults, 2: older adults), 2) start_time: start of
  follow-up, 3) end_time: end of follow-up, 4) time1: date for first
  serum collection, 5) time2: date for second serum collection, 6)
  time3: date for third serum collection, 7) HAI_titer_1: HAI titer for
  first serum collection, 8) HAI_titer_2: HAI titer for second serum
  collection, 9) HAI_titer_3: HAI titer for third serum collection.

- inputILI:

  The data for influenza activity used in the inference. The row number
  should match with the date in the inputdata.

- n_iteration:

  The number of iterations of the MCMC.

- burnin:

  The iteration for burn-in for MCMC.

- thinning:

  The number of thinning in MCMC.

- group_by:

  Optional formula specifying grouping variables (e.g., `~age_group`).
  When provided, independent MCMCs are fit for each combination of the
  grouping variables. The formula uses interaction semantics:
  `~age + vac` means all age-by-vac combinations. Returns a
  `seroreconstruct_multi` object.

- shared:

  Optional character vector specifying which parameters to share across
  groups when `group_by` is also provided. Measurement error parameters
  are always shared (they are a lab assay property, identical across
  groups). Valid values: `"error"` (measurement error only, the default
  when `shared` is non-NULL), `"boosting_waning"` (also share antibody
  boosting and waning across groups). When specified, a single joint
  MCMC is run with all groups pooled together, sharing the specified
  parameters while estimating group-specific infection probabilities.
  Returns a `seroreconstruct_joint` object.

- subject_ids:

  Optional vector (character, numeric, or factor) of subject
  identifiers, one per row of `inputdata`. When provided, stored in the
  fit object and used by
  [`plot_trajectory()`](https://timktsang.github.io/seroreconstruct/reference/plot_trajectory.md)
  to look up individuals by ID rather than row index. Example:
  `subject_ids = inputdata$household_id`.

## Value

A `seroreconstruct_fit` object (when `group_by` is `NULL`) or a
`seroreconstruct_multi` object (when `group_by` is provided). Use
[`summary()`](https://rdrr.io/r/base/summary.html) to extract model
estimates.

## Details

**Multi-season support:** If `inputdata` contains an optional integer
column named `season` (0-indexed, contiguous from 0 to `n_seasons - 1`),
the model fits season-specific infection risk and HAI protection
parameters. When no `season` column is present, all individuals are
assigned to a single season (`n_seasons = 1`) and behavior is identical
to previous versions. Validated with simulation recovery studies up to 7
seasons.

**Shared parameters:** When `shared` is provided together with
`group_by`, a single joint MCMC chain is run with all individuals
pooled. Measurement error and boosting/waning parameters are shared
across groups (informed by all data), while infection risk and HAI
protection parameters remain group-specific. This is more statistically
efficient than independent chains when groups share biological or
measurement properties.

**Single-group design:** When using `group_by` without `shared`,
independent MCMCs are fit for each group. To compare children vs adults,
fit each group separately using `group_by = ~age_group`.

**Current limitation:**
[`summary()`](https://rdrr.io/r/base/summary.html) is not yet
implemented for fits with `n_seasons > 1`. Multi-season posterior
samples are accessible directly from the fit object (e.g.,
`fit$posterior_model_parameter`).

## Examples

``` r
# \donttest{
a1 <- sero_reconstruct(inputdata, flu_activity,
                        n_iteration = 2000, burnin = 1000, thinning = 1)
#> 1000
#> MCMC complete in 28 seconds. Use summary() to view estimates.
summary(a1)
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
#>            2.38        1.44        3.65
#>            2.52        3.22        1.96
#>            7.34        4.61       10.71
#>            1.10        1.05        1.19
#>            5.50        4.43        6.78
#>            1.13        1.09        1.20
#>            0.21        0.18        0.24
#>            0.21        0.17        0.24
#>            0.16        0.11        0.21
#>            0.41        0.34        0.49
#>            0.28        0.23        0.34
#>            0.22        0.15        0.29
#>            1.02        0.81        1.28
#>            1.48        1.22        1.80
#>            0.78        0.53        1.05
#>            0.78        0.55        1.05
# }
```
