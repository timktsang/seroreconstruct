# Example of flu activity data

This is an example of the flu activity data used in the
`seroreconstruct` function. This data frame specifies the format of the
flu activity data.

## Usage

``` r
data(flu_activity)
```

## Format

A data frame with 1 variable, where each row represents a date, and it
should match the date in the input data:

- h1.activity:

  This is the influenza activity from surveillance data. It can be on a
  relative scale, as the model includes a scale parameter to estimate
  infection probability.

## See also

Other example_data:
[`para1`](https://timktsang.github.io/seroreconstruct/reference/para1.md),
[`para2`](https://timktsang.github.io/seroreconstruct/reference/para2.md)
