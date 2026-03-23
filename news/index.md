# Changelog

## seroreconstruct 1.1.0

### New features

- [`sero_reconstruct()`](https://timktsang.github.io/seroreconstruct/reference/sero_reconstruct.md)
  gains a `group_by` argument for subgroup analysis — fit independent
  MCMCs for age groups, vaccination status, or other strata defined by a
  one-sided formula (e.g., `~age_group`).

- [`sero_reconstruct()`](https://timktsang.github.io/seroreconstruct/reference/sero_reconstruct.md)
  gains a `shared` argument for joint models — run a single MCMC that
  shares measurement error (`"error"`) and/or antibody boosting/waning
  (`"boosting_waning"`) across groups while estimating group-specific
  infection risk. Measurement error is always shared when comparing
  groups of the same virus subtype.

- [`sero_reconstruct()`](https://timktsang.github.io/seroreconstruct/reference/sero_reconstruct.md)
  gains a `subject_ids` argument for ID-based individual lookup in
  [`plot_trajectory()`](https://timktsang.github.io/seroreconstruct/reference/plot_trajectory.md).

- **Multi-season support** — add a 0-indexed integer `season` column to
  input data; the model estimates season-specific infection risk and HAI
  protection parameters.

- **S3 classes** `seroreconstruct_fit`, `seroreconstruct_joint`, and
  `seroreconstruct_multi` with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods.

- **New plot functions:**
  [`plot_trajectory()`](https://timktsang.github.io/seroreconstruct/reference/plot_trajectory.md),
  [`plot_boosting()`](https://timktsang.github.io/seroreconstruct/reference/plot_boosting.md),
  [`plot_waning()`](https://timktsang.github.io/seroreconstruct/reference/plot_waning.md),
  [`plot_infection_prob()`](https://timktsang.github.io/seroreconstruct/reference/plot_infection_prob.md),
  [`plot_diagnostics()`](https://timktsang.github.io/seroreconstruct/reference/plot_diagnostics.md).

- **New table functions:**
  [`table_parameters()`](https://timktsang.github.io/seroreconstruct/reference/table_parameters.md),
  [`table_infections()`](https://timktsang.github.io/seroreconstruct/reference/table_infections.md).

- **[`simulate_data()`](https://timktsang.github.io/seroreconstruct/reference/simulate_data.md)**
  — generate synthetic HAI titer datasets for validation and power
  analysis.

### Deprecated

- [`output_model_estimate()`](https://timktsang.github.io/seroreconstruct/reference/output_model_estimate.md)
  is deprecated. Use
  [`table_parameters()`](https://timktsang.github.io/seroreconstruct/reference/table_parameters.md)
  and
  [`table_infections()`](https://timktsang.github.io/seroreconstruct/reference/table_infections.md)
  instead.

------------------------------------------------------------------------

## seroreconstruct 1.0.0

- Initial release. Core Bayesian MCMC framework for inferring influenza
  infection status, antibody dynamics, and individual infection risks
  from longitudinal HAI titer data.

- Based on Tsang TK et al. (2022) *Nat Commun* 13:1557.
  <https://doi.org/10.1038/s41467-022-29310-8>
