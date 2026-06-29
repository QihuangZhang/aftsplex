# Summarise a two-stage bootstrap SIMEX dose-response fit

Returns a tabular summary of the SIMEX-corrected dose-response curve and
its bootstrap confidence interval, together with the configuration that
produced it. The headline table reports the point estimate and 95%
percentile interval at quantiles of the exposure grid; the configuration
block reports the replicate yield, the exposure range and curve anchor,
and the spread of the Phase-1 `sigma_w_sq` across the replicates whose
calibration succeeded.

## Usage

``` r
# S3 method for class 'two_stage_bootstrap'
summary(object, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), ...)
```

## Arguments

- object:

  A
  [`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)
  fit.

- probs:

  Quantiles of `x_grid` at which to report the curve and its CI. Default
  `c(0.05, 0.25, 0.50, 0.75, 0.95)`.

- ...:

  Unused.

## Value

An object of class `"summary.two_stage_bootstrap"` with a `table`
element and summary fields. Has a `print` method.
