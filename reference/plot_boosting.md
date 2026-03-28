# Plot posterior boosting distributions

Draws violin plots of the posterior fold-rise in antibody titer after
infection, one violin per boosting parameter group. Matches the style of
Figure 1C in Tsang et al. (2022).

## Usage

``` r
plot_boosting(fit, cols = NULL, main = NULL, show_legend = TRUE, ...)
```

## Arguments

- fit:

  A `seroreconstruct_fit` or `seroreconstruct_joint` object. Only
  single-season fits are currently supported.

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
plot_boosting(fit)

# }
```
