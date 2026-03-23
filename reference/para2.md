# Example of parameter vector for the baseline HAI titer for the main model

This is an example of the parameter vector for the baseline HAI titer
for the main model used in the `seroreconstruct` function. This data
frame specifies the format of the parameter vector for the baseline HAI
titer for the main model.

## Usage

``` r
data(para2)
```

## Format

A numeric vector with 20 elements (for a single-season model, `S = 1`).
The general length is `20 * S` where `S` is the number of seasons.

- Elements 1–10:

  probability that the HAI titer is 0–9 for children

- Elements 11–20:

  probability that the HAI titer is 0–9 for adults

For multi-season models, this pattern repeats for each season.

## See also

Other example_data:
[`flu_activity`](https://timktsang.github.io/seroreconstruct/reference/flu_activity.md),
[`para1`](https://timktsang.github.io/seroreconstruct/reference/para1.md)
