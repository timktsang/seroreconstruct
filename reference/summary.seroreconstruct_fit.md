# Summary method for seroreconstruct_fit

Computes estimates of infection probabilities, boosting, waning, and
measurement error from a fitted MCMC object.

## Usage

``` r
# S3 method for class 'seroreconstruct_fit'
summary(object, period, ...)
```

## Arguments

- object:

  A `seroreconstruct_fit` object.

- period:

  Optional numeric vector of length 2 specifying the start and end of a
  season to compute infection probabilities. If omitted, the full
  follow-up period is used.

- ...:

  Additional arguments (ignored).

## Value

A `summary.seroreconstruct_fit` object with element `$table`.
