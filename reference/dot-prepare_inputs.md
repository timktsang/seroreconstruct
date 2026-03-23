# Prepare input data for MCMC or simulation

Reorders columns, adds padding columns, adjusts times for boosting
delay, and converts to matrix format expected by the C++ backend.

## Usage

``` r
.prepare_inputs(inputdata, inputILI)
```

## Arguments

- inputdata:

  Data frame with 9 required columns.

- inputILI:

  Data frame or matrix of influenza activity.

## Value

A list with prepared `inputdata` and `inputILI` matrices.
