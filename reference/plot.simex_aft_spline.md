# Plot a SIMEX-corrected AFT-spline fit

Convenience S3 method that draws the SIMEX-corrected curve via
[`plot_curves()`](https://qihuangzhang.github.io/aftsplex/reference/plot_curves.md).
By default the naive curve (`lambda = 0`) is overlaid as a
smoothing-bias diagnostic. Pass `truth` to add the dashed truth line,
and `ci` to add a percentile ribbon (typically from
[`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)).

## Usage

``` r
# S3 method for class 'simex_aft_spline'
plot(
  x,
  truth = NULL,
  ci = NULL,
  show_naive = TRUE,
  time_ratio = FALSE,
  title = "SIMEX-corrected AFT-spline dose-response",
  ...
)
```

## Arguments

- x:

  A
  [`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md)
  fit.

- truth:

  Optional numeric vector (length `length(x$x_grid)`) of the centred
  truth curve to overlay as a dashed black line.

- ci:

  Optional list with elements `lower` and `upper`.

- show_naive:

  Logical; overlay the naive (`lambda = 0`) curve? Default `TRUE`.

- time_ratio:

  Logical; plot the exponentiated curve (time ratio relative to the
  reference anchor) instead of the centred linear predictor? Default
  `FALSE`. Passed to
  [`plot_curves()`](https://qihuangzhang.github.io/aftsplex/reference/plot_curves.md).

- title:

  Plot title.

- ...:

  Passed through to
  [`plot_curves()`](https://qihuangzhang.github.io/aftsplex/reference/plot_curves.md).

## Value

A `ggplot` object.
