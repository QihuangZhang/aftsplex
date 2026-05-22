# Fit an AFT model with a natural-cubic-spline exposure

Thin wrapper around
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
that builds a formula of the form
`Surv(outcome_var, status_var) ~ ns(x_var, df = ...) + covariates` from
the column names supplied. Used three ways in the simulation: as the
oracle (on the true exposure), as the naive estimator (on the
GLS-combined surrogate without ME correction), and as the inner
workhorse of
[`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md).

## Usage

``` r
fit_aft_spline(
  data,
  x_var,
  covariates = character(0),
  df = 4,
  knots = NULL,
  outcome_var = "T_obs",
  status_var = "delta",
  dist = "lognormal",
  control = NULL
)
```

## Arguments

- data:

  A data frame containing the outcome, status, exposure, and covariate
  columns named below.

- x_var:

  Name of the exposure column.

- covariates:

  Character vector of additional covariate column names. Pass
  `character(0)` (the default) for an exposure-only model.

- df:

  Spline degrees of freedom (ignored if `knots` is supplied).

- knots:

  Optional numeric vector of interior knot locations.

- outcome_var:

  Name of the time-to-event column. Default `"T_obs"`.

- status_var:

  Name of the event indicator column. Default `"delta"`.

- dist:

  `survreg` distribution. Default `"lognormal"` matches the DGP in
  [`generate_aft_data()`](https://qihuangzhang.github.io/aftsplex/reference/generate_aft_data.md).

- control:

  Optional
  [`survival::survreg.control()`](https://rdrr.io/pkg/survival/man/survreg.control.html)
  object.

## Value

A `survreg` fitted-model object.

## Examples

``` r
sim <- generate_aft_data(n = 200, n_val = 100, seed = 1)
fit <- fit_aft_spline(sim$survival, x_var = "X_true",
                      covariates = c("V1","V2","V3","V4"))
coef(fit)
#>         (Intercept) ns(X_true, df = 4)1 ns(X_true, df = 4)2 ns(X_true, df = 4)3 
#>          2.95825639          1.14551074          0.79396526          1.45436664 
#> ns(X_true, df = 4)4                  V1                  V2                  V3 
#>          0.38171431          0.11330179         -0.05156747          0.50880591 
#>                  V4 
#>         -0.11017939 
```
