# Validate simulation parameters

Validate simulation parameters

## Usage

``` r
.validate_simulation_params(para1, para2, n_seasons = 1L, n_groups = 3L)
```

## Arguments

- para1:

  Numeric vector of active model parameters (length
  `6 + 4 * n_seasons`).

- para2:

  Numeric vector of baseline HAI titer distribution (length
  `20 * n_seasons`).

- n_seasons:

  Number of seasons. Default 1.
