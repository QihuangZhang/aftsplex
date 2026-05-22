# GLS-combine surrogate columns into a single calibrated exposure

Applies the linear back-transform `W_tilde = (W - alpha_0) / alpha_1`
and the GLS combiner `omega = Sigma_tilde^-1 1 / (1' Sigma_tilde^-1 1)`
to produce a single calibrated exposure with conditional variance
`sigma_w_sq = (1' Sigma_tilde^-1 1)^-1`.

## Usage

``` r
gls_combine(W, calibration)
```

## Arguments

- W:

  Numeric matrix of surrogate columns, columns ordered as in
  `calibration$W_cols`.

- calibration:

  Output of
  [`fit_me_calibration()`](https://qihuangzhang.github.io/aftsplex/reference/fit_me_calibration.md).

## Value

A list with elements

- `W_bar` numeric vector of length `nrow(W)`,

- `sigma_w_sq` scalar conditional variance,

- `omega` numeric vector of GLS weights (length `J`),

- `Sigma_tilde` `J x J` back-transformed error covariance.

## Examples

``` r
sim <- generate_aft_data(n = 100, n_val = 200, seed = 1)
cal <- fit_me_calibration(sim$validation)
W   <- as.matrix(sim$survival[, cal$W_cols])
out <- gls_combine(W, cal)
head(out$W_bar)
#> [1]  9.435341 11.776481  9.556977 12.020854 11.680636  9.901770
```
