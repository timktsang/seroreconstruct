# Simulation of the dataset of the Bayesian model

The function to simulate the dataset, for validation or other purpose.

## Usage

``` r
simulate_data(inputdata, inputILI, para1, para2, n_groups = 3L)
```

## Arguments

- inputdata:

  The data with the same format that for running MCMC, in dataframe
  format.

- inputILI:

  The data for influenza activity used in the inference. The row number
  should match with the date in the inputdata.

- para1:

  Numeric vector of active model parameters. Length depends on the
  number of seasons `S` (determined by the `season` column in
  `inputdata`, default `S = 1`):

  - Elements 1–6 (shared): 1) random measurement error, 2) 2-fold
    error, 3) boosting for children (log2), 4) waning for children
    (log2), 5) boosting for adults (log2), 6) waning for adults (log2).

  - Elements 7 to `6 + 3*S` (per-season): infection risk scale
    parameters for children, adults, and older adults, repeated for each
    season.

  - Elements `6 + 3*S + 1` to `6 + 4*S` (per-season): log risk ratio of
    2-fold increase in baseline HAI titer, one per season.

  Total length: `6 + (G + 1)*S` where `G` is `n_groups` (e.g., 10 for
  G=3 S=1, 34 for G=3 S=7). See
  [`para1`](https://timktsang.github.io/seroreconstruct/reference/para1.md)
  for an example with `G = 3, S = 1`.

- para2:

  Numeric vector for baseline HAI titer distributions. Length `20 * S`:
  for each season, 10 probabilities for children (HAI titer levels 0–9)
  followed by 10 probabilities for adults. See
  [`para2`](https://timktsang.github.io/seroreconstruct/reference/para2.md)
  for an example with `S = 1`.

- n_groups:

  Number of groups for infection risk parameters (default 3 for the
  standard 3-age-group model).

## Value

A simulated data based on the input parameter vectors, with the format
equal to the input data.

## Examples

``` r
simulated <- simulate_data(inputdata, flu_activity, para1, para2)
```
