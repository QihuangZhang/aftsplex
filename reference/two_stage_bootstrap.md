# Two-stage nonparametric bootstrap for the SIMEX-corrected curve

For each replicate `r`:

1.  resample validation participants;

2.  refit the multivariate ME calibration via
    [`fit_me_calibration()`](https://qihuangzhang.github.io/aftsplex/reference/fit_me_calibration.md);

3.  resample main-study participants;

4.  form the bootstrap-combined surrogate via
    [`gls_combine()`](https://qihuangzhang.github.io/aftsplex/reference/gls_combine.md)
    under the refit calibration;

5.  run
    [`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md)
    with the refit `sigma_w_sq` and store the corrected curve.

## Usage

``` r
two_stage_bootstrap(
  survival,
  validation,
  x_var = "X_true",
  covariates = character(0),
  v_ref = NULL,
  surrogate_pattern = "^W[0-9]+$",
  df = 4,
  knots = NULL,
  outcome_var = "T_obs",
  status_var = "delta",
  lambda = c(0.5, 1, 1.5, 2),
  B = 50,
  R = 200,
  x_grid = NULL,
  dist = "lognormal",
  verbose = FALSE
)
```

## Arguments

- survival:

  Main-study data frame; must contain the columns named in `x_var` (if a
  single calibrated exposure already exists), or the surrogate columns
  matched by `surrogate_pattern`, plus the outcome, status, and
  covariate columns.

- validation:

  Validation data frame containing the truth column named in `x_var` and
  the surrogate columns matched by `surrogate_pattern`.

- x_var:

  Name of the exposure column. For simulation use the truth column name
  (`"X_true"`); for the standard workflow this is the same name used in
  both `validation` and the bootstrap-recombined `survival` data.

- covariates:

  Character vector of confounder column names.

- v_ref:

  Named numeric vector of reference covariate values matching
  `covariates`. `NULL` if no covariates.

- surrogate_pattern:

  Regex used to identify surrogate columns; default `"^W[0-9]+$"`.

- df, knots, lambda, B, dist:

  Passed to
  [`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md).

- outcome_var, status_var:

  Passed to
  [`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md).

- R:

  Number of outer bootstrap replicates.

- x_grid:

  Numeric exposure grid for the curve. If `NULL`, taken as 100 equally
  spaced points between the 2nd and 98th percentile of the full-sample
  combined surrogate.

- verbose:

  Print progress bar?

## Value

A list with elements

- `x_grid` exposure grid,

- `f_hat` single-fit SIMEX point estimate (full-sample, no resampling),

- `bs_mean` mean of the bootstrap curves,

- `curves` matrix of bootstrap curves (`length(x_grid)` rows by `R`
  cols),

- `median, lower, upper` percentile-CI summaries (2.5%, 50%, 97.5%),

- `lower_basic, upper_basic` reverse-percentile (basic) CI bounds,

- `f_bc` bias-corrected point estimate `2*f_hat - bs_mean`,

- `sigma_w_sq_boot, omega_boot` per-replicate Phase-1 summaries,

- `R_effective` number of replicates that produced a curve.

Empirically (see
[`vignette("coverage-diagnostics")`](https://qihuangzhang.github.io/aftsplex/articles/coverage-diagnostics.md))
the percentile CI is the primary interval; the basic CI is kept for
completeness but does not improve coverage on simulated
diminishing-returns curves because the dominant residual bias is from
spline smoothing rather than SIMEX.

## Details

Pointwise quantiles across `r` give CIs that propagate both Phase-1 and
Phase-2 uncertainty - this addresses the plug-in-variance critique that
the conditional `sigma_w_sq` alone yields anti-conservative intervals
when the validation sample is small relative to the main study.

## Examples

``` r
# \donttest{
sim <- generate_aft_data(n = 500, n_val = 200, seed = 1)
boot <- two_stage_bootstrap(
  survival   = sim$survival,
  validation = sim$validation,
  x_var      = "X_true",
  covariates = c("V1","V2","V3","V4"),
  v_ref      = c(V1=30, V2=30, V3=0, V4=0),
  lambda     = c(0.5, 1, 1.5, 2), B = 5, R = 10
)
head(boot$f_hat)
#>          1          2          3          4          5          6 
#> 0.00000000 0.01912404 0.03855214 0.05830025 0.07838434 0.09882034 
# }
```
