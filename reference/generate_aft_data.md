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
#> 1  156.71416     1  8.599207  9.042557  9.435666  8.025849 32.04701 35.37220  1
#> 2  110.70489     1 10.410639  9.409792 10.765985 10.894197 38.44437 39.47827  0
#> 3  333.80655     0  8.131478  8.497433  9.839325  8.418137 37.93294 26.98501  0
#> 4  443.05466     1 13.567156 14.752046 15.172707 13.659899 28.34546 28.04566  1
#> 5   73.07738     1 10.736802  9.950293 10.977064  8.352338 18.57382 27.91889  0
#> 6 1150.60870     1  8.165377  8.262311  9.389844  6.764756 42.48831 28.12171  0
#>   V4
#> 1  1
#> 2  1
#> 3  0
#> 4  1
#> 5  1
#> 6  0
```
