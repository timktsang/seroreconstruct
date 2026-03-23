# Example of parameter vector for the main model

This is an example of the parameter vector for the main model used in
the `seroreconstruct` function. This data frame specifies the format of
the parameter vector for the main model.

## Usage

``` r
data(para1)
```

## Format

A numeric vector with 10 elements (for a single-season model, `S = 1`).
The general length is `6 + 4*S` where `S` is the number of seasons.

- Elements 1–6 (shared):

  1\) random measurement error, 2) 2-fold error, 3) boosting for
  children (log2), 4) waning for children (log2), 5) boosting for adults
  (log2), 6) waning for adults (log2).

- Elements 7–9 (per-season):

  infection risk scale parameters for children, adults, and older adults
  (3 per season).

- Element 10 (per-season):

  log risk ratio of 2-fold increase in baseline HAI titer (1 per
  season).

## See also

Other example_data:
[`flu_activity`](https://timktsang.github.io/seroreconstruct/reference/flu_activity.md),
[`para2`](https://timktsang.github.io/seroreconstruct/reference/para2.md)
