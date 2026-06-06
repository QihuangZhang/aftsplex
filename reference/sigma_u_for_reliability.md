# Build a measurement-error covariance targeting a given reliability

Returns a `Sigma_u` (via
[`build_sigma_u()`](https://qihuangzhang.github.io/aftsplex/reference/build_sigma_u.md))
whose surrogate-measurement error variances are scaled so the **mean
per-surrogate reliability** equals `reliability`. Surrogate reliability
is the intraclass-style ratio `var_x / (var_x + var_U_j)`; the relative
surrogate qualities in `ratios` are preserved (so the surrogates stay
unequal for the per-surrogate ladder), and a single multiplicative scale
is solved by [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html)
to hit the target average. Use it to sweep measurement quality in a
simulation, e.g. `reliability` in `c(0.5, 0.65, 0.8)`.

## Usage

``` r
sigma_u_for_reliability(reliability, var_x, ratios = c(2, 2.5, 3), rho = 0.4)
```

## Arguments

- reliability:

  Target mean per-surrogate reliability in `(0, 1)`.

- var_x:

  Variance of the latent exposure `X` (e.g. `x_sd^2` passed to
  [`generate_aft_data()`](https://qihuangzhang.github.io/aftsplex/reference/generate_aft_data.md)).

- ratios:

  Relative marginal error variances across the surrogates; their scale
  is solved for, only their ratios matter. Default `c(2.0, 2.5, 3.0)`
  matches
  [`build_sigma_u()`](https://qihuangzhang.github.io/aftsplex/reference/build_sigma_u.md).

- rho:

  Common pairwise error correlation, passed to
  [`build_sigma_u()`](https://qihuangzhang.github.io/aftsplex/reference/build_sigma_u.md).

## Value

A `J x J` covariance matrix suitable for the `Sigma_u` argument of
[`generate_aft_data()`](https://qihuangzhang.github.io/aftsplex/reference/generate_aft_data.md).

## Examples

``` r
S <- sigma_u_for_reliability(0.65, var_x = 75^2)
# mean per-surrogate reliability is ~0.65:
mean(75^2 / (75^2 + diag(S)))
#> [1] 0.65
```
