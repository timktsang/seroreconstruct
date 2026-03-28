# Plot posterior waning curves

Shows the fraction of peak antibody remaining over time since infection,
with posterior median and 95% credible band for each waning parameter
group. Matches the style of Figure 1D in Tsang et al. (2022).

## Usage

``` r
plot_waning(fit, days = 400, cols = NULL, main = NULL, show_legend = TRUE, ...)
```

## Arguments

- fit:

  A `seroreconstruct_fit` or `seroreconstruct_joint` object. Only
  single-season fits are currently supported.

- days:

  Maximum number of days to plot on the x-axis. Default 400.

- cols:

  Optional character vector of colors, one per group.

- main:

  Optional plot title.

- show_legend:

  Logical; whether to draw a legend. Default `TRUE`.

- ...:

  Additional graphical parameters passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisible `NULL`. Called for its side effect of producing a plot.

## Examples

``` r
# \donttest{
fit <- sero_reconstruct(inputdata, flu_activity,
                        n_iteration = 2000, burnin = 1000, thinning = 1)
#> 1000
#> MCMC complete in 27 seconds. Use summary() to view estimates.
plot_waning(fit)

# }
```
