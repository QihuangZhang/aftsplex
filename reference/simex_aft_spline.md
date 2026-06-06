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
  x_ref = NULL,
  support_probs = NULL,
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

- x_ref:

  Optional scalar reference exposure at which to anchor the centred
  curves (so the dose-response is read as relative to `x_ref`, e.g. a
  clinically meaningful value). `NULL` (default) anchors at the lower
  grid boundary `x_grid[1]`, preserving the original behaviour.

- support_probs:

  Optional length-2 numeric of probabilities defining the trustworthy
  exposure support as quantiles of `data[[x_var]]` (e.g.
  `c(0.05, 0.95)`). When set, grid points outside the support are
  flagged in the returned `in_support` vector and a one-shot
  [`warning()`](https://rdrr.io/r/base/warning.html) notes that they are
  spline extrapolations. `NULL` (default) disables the check.

- verbose:

  Print progress bar and non-convergence count?

- maxiter:

  `survreg` iteration cap (passed via
  [`survival::survreg.control()`](https://rdrr.io/pkg/survival/man/survreg.control.html)).

## Value

An object of class `"simex_aft_spline"`: a list with elements

- `x_grid` the exposure grid,

- `curve_simex` the centred SIMEX-corrected linear predictor,

- `curve_naive` the centred naive linear predictor (lambda = 0),

- `curves_lambda` matrix of centred lambda-fits (columns = lambda,
  including 0),

- `lambda` the lambda values used (with 0 in column 1),

- `in_support` logical vector aligned to `x_grid` (or `NULL` if
  `support_probs` was not supplied),

- `x_ref` the reference anchor used (or `NULL`),

- `call`, `n`, `df`, `dist`, `covariates`, `B`, `sigma_w_sq`, `n_fail`
  configuration and diagnostics used by
  [`summary.simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/summary.simex_aft_spline.md)
  and
  [`plot.simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/plot.simex_aft_spline.md).

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
summary(s)
#> SIMEX-corrected AFT-spline dose-response
#> 
#> Call:
#>   simex_aft_spline(data = dat, x_var = "W_bar", sigma_w_sq = g$sigma_w_sq, 
#>     covariates = c("V1", "V2", "V3", "V4"), v_ref = c(V1 = 30, 
#>         V2 = 30, V3 = 0, V4 = 0), lambda = c(0.5, 1, 1.5, 2), 
#>     B = 10)
#> 
#> Configuration:
#>   n observations:  500                Spline df:    4
#>   Covariates:      V1, V2, V3, V4     Distribution: lognormal
#>   Lambda grid:     0.0, 0.5, 1.0, 1.5, 2.0
#>   Inner B:         10                 sigma_w_sq:   1.2659
#> 
#> Centred dose-response (anchored at min(x_grid)):
#>               5%   25%    50%    75%    95%
#> x          3.617 6.666 10.478 14.289 17.339
#> Naive      0.095 0.514  1.007  0.960  1.621
#> SIMEX      0.111 0.592  1.199  1.224  2.208
#> Correction 0.016 0.078  0.192  0.264  0.587
# }
```
