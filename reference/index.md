# Package index

## Model Fitting

Core functions for running the MCMC and extracting results

- [`sero_reconstruct()`](https://timktsang.github.io/seroreconstruct/reference/sero_reconstruct.md)
  : Run the MCMC for the Bayesian model
- [`summary(`*`<seroreconstruct_fit>`*`)`](https://timktsang.github.io/seroreconstruct/reference/summary.seroreconstruct_fit.md)
  : Summary method for seroreconstruct_fit
- [`summary(`*`<seroreconstruct_joint>`*`)`](https://timktsang.github.io/seroreconstruct/reference/summary.seroreconstruct_joint.md)
  : Summary method for seroreconstruct_joint
- [`summary(`*`<seroreconstruct_multi>`*`)`](https://timktsang.github.io/seroreconstruct/reference/summary.seroreconstruct_multi.md)
  : Summary method for seroreconstruct_multi

## Visualization

Publication-ready plots

- [`plot_diagnostics()`](https://timktsang.github.io/seroreconstruct/reference/plot_diagnostics.md)
  : MCMC diagnostic plots
- [`plot_trajectory()`](https://timktsang.github.io/seroreconstruct/reference/plot_trajectory.md)
  : Plot antibody trajectory with model fit
- [`plot_boosting()`](https://timktsang.github.io/seroreconstruct/reference/plot_boosting.md)
  : Plot posterior boosting distributions
- [`plot_waning()`](https://timktsang.github.io/seroreconstruct/reference/plot_waning.md)
  : Plot posterior waning curves
- [`plot_infection_prob()`](https://timktsang.github.io/seroreconstruct/reference/plot_infection_prob.md)
  : Plot infection probabilities (forest plot)

## Tables

Summary tables for model output

- [`table_parameters()`](https://timktsang.github.io/seroreconstruct/reference/table_parameters.md)
  : Summary table of model parameters with credible intervals
- [`table_infections()`](https://timktsang.github.io/seroreconstruct/reference/table_infections.md)
  : Per-individual infection estimates

## Simulation

Generate synthetic data

- [`simulate_data()`](https://timktsang.github.io/seroreconstruct/reference/simulate_data.md)
  : Simulation of the dataset of the Bayesian model

## Data

Bundled example datasets

- [`inputdata`](https://timktsang.github.io/seroreconstruct/reference/inputdata.md)
  : Example of input data
- [`flu_activity`](https://timktsang.github.io/seroreconstruct/reference/flu_activity.md)
  : Example of flu activity data
- [`para1`](https://timktsang.github.io/seroreconstruct/reference/para1.md)
  : Example of parameter vector for the main model
- [`para2`](https://timktsang.github.io/seroreconstruct/reference/para2.md)
  : Example of parameter vector for the baseline HAI titer for the main
  model

## S3 Methods

Print and subset methods

- [`print(`*`<seroreconstruct_fit>`*`)`](https://timktsang.github.io/seroreconstruct/reference/print.seroreconstruct_fit.md)
  : Print method for seroreconstruct_fit
- [`print(`*`<seroreconstruct_joint>`*`)`](https://timktsang.github.io/seroreconstruct/reference/print.seroreconstruct_joint.md)
  : Print method for seroreconstruct_joint
- [`print(`*`<seroreconstruct_multi>`*`)`](https://timktsang.github.io/seroreconstruct/reference/print.seroreconstruct_multi.md)
  : Print method for seroreconstruct_multi
- [`print(`*`<summary.seroreconstruct_fit>`*`)`](https://timktsang.github.io/seroreconstruct/reference/print.summary.seroreconstruct_fit.md)
  : Print method for summary.seroreconstruct_fit
- [`print(`*`<summary.seroreconstruct_joint>`*`)`](https://timktsang.github.io/seroreconstruct/reference/print.summary.seroreconstruct_joint.md)
  : Print method for summary.seroreconstruct_joint
- [`print(`*`<summary.seroreconstruct_multi>`*`)`](https://timktsang.github.io/seroreconstruct/reference/print.summary.seroreconstruct_multi.md)
  : Print method for summary.seroreconstruct_multi
- [`` `[[`( ``*`<seroreconstruct_multi>`*`)`](https://timktsang.github.io/seroreconstruct/reference/sub-sub-.seroreconstruct_multi.md)
  : Subset a seroreconstruct_multi object

## Deprecated

- [`output_model_estimate()`](https://timktsang.github.io/seroreconstruct/reference/output_model_estimate.md)
  : Extract the model estimates from the fitted MCMC
