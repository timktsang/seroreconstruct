# Summary method for seroreconstruct_multi

Computes estimates for each group and combines into a single table.

## Usage

``` r
# S3 method for class 'seroreconstruct_multi'
summary(object, period, ...)
```

## Arguments

- object:

  A `seroreconstruct_multi` object.

- period:

  Optional numeric vector of length 2 specifying the start and end of a
  season to compute infection probabilities.

- ...:

  Additional arguments (ignored).

## Value

A `summary.seroreconstruct_multi` object with element `$table`.
