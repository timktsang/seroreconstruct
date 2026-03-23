# seroreconstruct 1.1.0

## New features

* `sero_reconstruct()` gains a `group_by` argument for subgroup analysis — fit
  independent MCMCs for age groups, vaccination status, or other strata defined
  by a one-sided formula (e.g., `~age_group`).

* `sero_reconstruct()` gains a `shared` argument for joint models — run a single
  MCMC that shares measurement error (`"error"`) and/or antibody
  boosting/waning (`"boosting_waning"`) across groups while estimating
  group-specific infection risk. Measurement error is always shared when
  comparing groups of the same virus subtype.

* `sero_reconstruct()` gains a `subject_ids` argument for ID-based individual
  lookup in `plot_trajectory()`.

* **Multi-season support** — add a 0-indexed integer `season` column to input
  data; the model estimates season-specific infection risk and HAI protection
  parameters.

* **S3 classes** `seroreconstruct_fit`, `seroreconstruct_joint`, and
  `seroreconstruct_multi` with `print()` and `summary()` methods.

* **New plot functions:** `plot_trajectory()`, `plot_boosting()`,
  `plot_waning()`, `plot_infection_prob()`, `plot_diagnostics()`.

* **New table functions:** `table_parameters()`, `table_infections()`.

* **`simulate_data()`** — generate synthetic HAI titer datasets for validation
  and power analysis.

## Deprecated

* `output_model_estimate()` is deprecated. Use `table_parameters()` and
  `table_infections()` instead.

---

# seroreconstruct 1.0.0

* Initial release. Core Bayesian MCMC framework for inferring influenza
  infection status, antibody dynamics, and individual infection risks from
  longitudinal HAI titer data.

* Based on Tsang TK et al. (2022) *Nat Commun* 13:1557.
  <https://doi.org/10.1038/s41467-022-29310-8>
