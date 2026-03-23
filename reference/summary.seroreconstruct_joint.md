# Summary method for seroreconstruct_joint

Computes shared parameter estimates and per-group infection
probabilities from a joint fit with shared parameters.

## Usage

``` r
# S3 method for class 'seroreconstruct_joint'
summary(object, period, ...)
```

## Arguments

- object:

  A `seroreconstruct_joint` object.

- period:

  Optional numeric vector of length 2 specifying the start and end of a
  season to compute infection probabilities.

- ...:

  Additional arguments (ignored).

## Value

A `summary.seroreconstruct_joint` object with element `$table`.
