# Summarise a SIMEX-corrected AFT-spline fit

Returns a tabular summary of the SIMEX point estimate together with the
configuration that produced it. The headline table reports the centred
Naive and SIMEX dose-response (and their difference, the SIMEX
correction) at quantiles of the exposure grid.

## Usage

``` r
# S3 method for class 'simex_aft_spline'
summary(object, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), ...)
```

## Arguments

- object:

  A
  [`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md)
  fit.

- probs:

  Quantiles of `x_grid` at which to report the centred curves. Default
  `c(0.05, 0.25, 0.50, 0.75, 0.95)`.

- ...:

  Unused.

## Value

An object of class `"summary.simex_aft_spline"` with a `table` element
and the configuration fields. Has a `print` method.
