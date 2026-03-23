# MCMC diagnostic plots

Produces trace plots and posterior density plots for each model
parameter. Trace plots show the MCMC chain with the posterior mean (red
dashed line). Density plots show the marginal posterior with 95%
credible interval bounds (blue dashed lines).

## Usage

``` r
plot_diagnostics(fit, params = NULL)
```

## Arguments

- fit:

  A `seroreconstruct_fit`, `seroreconstruct_joint`, or
  `seroreconstruct_multi` object.

- params:

  Optional character vector of parameter names to plot. If `NULL`
  (default), all parameters are plotted. Use
  `table_parameters(fit)$Parameter` to see available names.

## Value

Invisible `NULL`. Called for its side effect of producing plots.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- sero_reconstruct(inputdata, flu_activity,
                        n_iteration = 2000, burnin = 1000, thinning = 1)
plot_diagnostics(fit)
} # }
```
