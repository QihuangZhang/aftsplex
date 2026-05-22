# Predict the centred dose-response curve at a reference covariate setting

Returns the linear predictor of a fitted
[`fit_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/fit_aft_spline.md)
model on an exposure grid, holding the covariates fixed at `v_ref`.
Subtracting the first element produces the centred curve used throughout
the manuscript.

## Usage

``` r
predict_curve(fit, x_grid, v_ref = NULL, x_var)
```

## Arguments

- fit:

  A `survreg` object returned by
  [`fit_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/fit_aft_spline.md).

- x_grid:

  Numeric vector of exposure values to predict at.

- v_ref:

  Named numeric vector of reference covariate values, with names
  matching the `covariates` used when fitting. Pass `NULL` if the model
  has no covariates.

- x_var:

  Name of the exposure column (must match `fit_aft_spline`).

## Value

A numeric vector of length `length(x_grid)` (the uncentred linear
predictor; subtract `out[1]` to centre at the lower grid boundary).

## Examples

``` r
sim <- generate_aft_data(n = 300, n_val = 100, seed = 1)
fit <- fit_aft_spline(sim$survival, x_var = "X_true",
                      covariates = c("V1","V2","V3","V4"))
x_grid <- seq(7, 13, length.out = 50)
lp <- predict_curve(fit, x_grid,
                    v_ref = c(V1=30, V2=30, V3=0, V4=0),
                    x_var = "X_true")
```
