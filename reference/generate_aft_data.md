# Generate one Monte Carlo replicate of the AFT-spline-SIMEX scenario

Produces a main-study survival sample and an external validation sample
from the same data-generating process described in the manuscript
appendix. The latent log-exposure is normal, the dose-response is the
saturating curve
[`f_true()`](https://qihuangzhang.github.io/aftsplex/reference/f_true.md),
the survival times follow a log-normal AFT model with four confounders
(V1-V4), and the surrogate measurements are multivariate-normal around
the latent truth with covariance `Sigma_u`.

## Usage

``` r
generate_aft_data(
  n = 2000,
  n_val = 500,
  Sigma_u = build_sigma_u(),
  mu = 4,
  sigma_T = 0.8,
  gamma = c(0.1, -0.05, 0.3, -0.2),
  cens_rate = 5e-04,
  seed = NULL
)
```

## Arguments

- n:

  Main-study sample size.

- n_val:

  Validation sample size.

- Sigma_u:

  `J x J` measurement-error covariance (default 3x3 with marginal
  variances `(2.0, 2.5, 3.0)` and `rho = 0.4`).

- mu:

  AFT intercept.

- sigma_T:

  AFT log-scale standard deviation.

- gamma:

  Numeric vector of length 4: confounder coefficients in the linear
  predictor (`V1, V2, V3, V4`).

- cens_rate:

  Rate of the exponential censoring distribution. Smaller values give
  heavier event experience. The default `0.0005` yields ~76% events
  under the other defaults.

- seed:

  Optional integer seed.

## Value

A list with elements `survival` (main-study data frame), `validation`
(validation data frame with paired truth + surrogates), and `truth` (a
list of the parameter values used).

## Details

The returned `survival` data frame is in the shape expected by
[`fit_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/fit_aft_spline.md)
and
[`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)
with their defaults (`outcome_var = "T_obs"`, `status_var = "delta"`,
`covariates = c("V1", "V2", "V3", "V4")`, `surrogate_cols = NULL` -\>
auto-detect via
[`extract_surrogate_cols()`](https://qihuangzhang.github.io/aftsplex/reference/extract_surrogate_cols.md)).

## Examples

``` r
sim <- generate_aft_data(n = 200, n_val = 100, seed = 1)
names(sim)
#> [1] "survival"   "validation" "truth"     
head(sim$survival)
#>        T_obs delta    X_true        W1        W2        W3       V1       V2 V3
#> 1  156.71416     1  8.599207  8.357699  7.655838  9.125222 32.04701 35.37220  1
#> 2  110.70489     1 10.410639  9.461024 11.088390 10.384565 38.44437 39.47827  0
#> 3  333.80655     0  8.131478  6.858229  6.904194  8.057625 37.93294 26.98501  0
#> 4  443.05466     1 13.567156 12.850575 11.713557 13.364571 28.34546 28.04566  1
#> 5   73.07738     1 10.736802 10.992466 10.777710 13.245776 18.57382 27.91889  0
#> 6 1150.60870     1  8.165377  7.734338  7.117875  9.644362 42.48831 28.12171  0
#>   V4
#> 1  1
#> 2  1
#> 3  0
#> 4  1
#> 5  1
#> 6  0
```
