# SIMEX-correct the AFT-spline dose-response curve

Implements simulation-extrapolation (Cook and Stefanski 1994; Carroll,
Kuchenhoff, Lombard and Stefanski 1996) for a natural-cubic-spline AFT
model fitted by
[`fit_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/fit_aft_spline.md).
For each `lambda` in the supplied grid, draws `B` perturbed surrogates
with additional noise `sqrt(lambda * sigma_w_sq) * N(0, 1)`, refits,
averages the centred linear predictor across `B`, then extrapolates to
`lambda = -1` via a pointwise quadratic OLS.

## Usage

``` r
simex_aft_spline(
  data,
  x_var,
  sigma_w_sq,
  covariates = character(0),
  v_ref = NULL,
  df = 4,
  knots = NULL,
  outcome_var = "T_obs",
  status_var = "delta",
  lambda = c(0.5, 1, 1.5, 2),
  B = 50,
  dist = "lognormal",
  x_grid = NULL,
  verbose = FALSE,
  maxiter = 200
)
```

## Arguments

- data:

  A data frame containing the columns named in `x_var`, `covariates`,
  `outcome_var`, and `status_var`.

- x_var:

  Name of the exposure column.

- sigma_w_sq:

  Scalar conditional variance of the calibrated exposure (typically
  `gls_combine(...)$sigma_w_sq`).

- covariates:

  Character vector of confounder column names.

- v_ref:

  Named numeric vector of reference covariate values matching
  `covariates`. `NULL` if no covariates.

- df, knots, dist:

  Passed to
  [`fit_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/fit_aft_spline.md).

- outcome_var, status_var:

  Passed to
  [`fit_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/fit_aft_spline.md).

- lambda:

  Numeric grid of SIMEX lambda values (excluding 0; 0 is appended
  internally as the naive fit). Default `c(0.5, 1, 1.5, 2)`.

- B:

  Number of inner replicates per lambda.

- x_grid:

  Numeric grid of exposure values at which to evaluate the corrected
  curve. Defaults to 100 equally spaced points spanning `data[[x_var]]`.

- verbose:

  Print progress bar and non-convergence count?

- maxiter:

  `survreg` iteration cap (passed via
  [`survival::survreg.control()`](https://rdrr.io/pkg/survival/man/survreg.control.html)).

## Value

A list with elements

- `x_grid` the exposure grid,

- `curve_simex` the centred SIMEX-corrected linear predictor,

- `curves_lambda` matrix of centred lambda-fits (columns = lambda,
  including 0),

- `lambda` the lambda values used (with 0 in column 1).

## Details

Inner `survreg` fits are wrapped in convergence handling: fits that fail
to converge within `maxiter` iterations are dropped from the average via
`na.rm = TRUE`, and the count is reported if `verbose = TRUE`.

## Examples

``` r
# \donttest{
sim <- generate_aft_data(n = 500, n_val = 200, seed = 1)
cal <- fit_me_calibration(sim$validation)
W   <- as.matrix(sim$survival[, cal$W_cols])
g   <- gls_combine(W, cal)
dat <- sim$survival
dat$W_bar <- g$W_bar
s <- simex_aft_spline(dat, x_var = "W_bar",
                      sigma_w_sq = g$sigma_w_sq,
                      covariates = c("V1","V2","V3","V4"),
                      v_ref = c(V1=30, V2=30, V3=0, V4=0),
                      lambda = c(0.5, 1, 1.5, 2), B = 10)
head(s$curve_simex)
#>          1          2          3          4          5          6 
#> 0.00000000 0.02254499 0.04499498 0.06737831 0.08972246 0.11205485 
# }
```
