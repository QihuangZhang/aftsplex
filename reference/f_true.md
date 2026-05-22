# True dose-response used by the simulation DGP

Saturating exponential `beta * (1 - exp(-alpha * pmax(x - x_min, 0)))`.
This is the function the simulation tries to recover; it serves as the
ground truth in
[`generate_aft_data()`](https://qihuangzhang.github.io/aftsplex/reference/generate_aft_data.md)
and is exported so users can overlay it on plots.

## Usage

``` r
f_true(x, beta = 0.6, alpha = 0.4, x_min = 6)
```

## Arguments

- x:

  Numeric vector of exposure values.

- beta:

  Saturation amplitude (asymptote as `x` grows large).

- alpha:

  Curvature parameter (larger = steeper initial rise).

- x_min:

  Threshold below which the response is flat at zero.

## Value

Numeric vector the same length as `x`.

## Examples

``` r
x <- seq(0, 20, length.out = 100)
plot(x, f_true(x), type = "l")

```
