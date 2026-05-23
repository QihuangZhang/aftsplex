# Overlay fitted dose-response curves with the truth and an optional CI

Builds a faceted ggplot2 figure comparing a named list of fitted curves
against the truth, with an optional bootstrap confidence ribbon. Uses
the
[`method_colors()`](https://qihuangzhang.github.io/aftsplex/reference/method_colors.md)
palette by default.

## Usage

``` r
plot_curves(
  x_grid,
  truth,
  fits,
  ci = NULL,
  title = "AFT + spline dose-response",
  x_lab = "Exposure",
  y_lab = "Centred linear predictor (log time ratio)",
  colors = method_colors()
)
```

## Arguments

- x_grid:

  Numeric exposure grid.

- truth:

  Numeric vector of truth values (already centred at the first grid
  point), the same length as `x_grid`. Pass `NULL` to omit the truth
  line (e.g. on real data where the truth is unknown).

- fits:

  Named list of numeric vectors of length `length(x_grid)`.

- ci:

  Optional list with elements `lower` and `upper`, each a numeric vector
  of length `length(x_grid)`.

- title:

  Plot title.

- x_lab, y_lab:

  Axis labels.

- colors:

  Optional named character vector of colours; defaults to
  [`method_colors()`](https://qihuangzhang.github.io/aftsplex/reference/method_colors.md)
  (matching `Naive`, `Oracle`, `SIMEX`).

## Value

A `ggplot` object.

## Details

Requires the suggested package `ggplot2` to be installed.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  x <- seq(6, 14, length.out = 50)
  truth <- f_true(x) - f_true(x[1])
  plot_curves(x, truth, list(Oracle = truth))
}

```
