# Semantic method-to-colour mapping for Naive / Oracle / SIMEX comparisons

Maps the three estimator names used throughout the package to colours
drawn from
[`stagill_palette()`](https://qihuangzhang.github.io/aftsplex/reference/stagill_palette.md).
Pass to
[`ggplot2::scale_color_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Usage

``` r
method_colors()
```

## Value

A named character vector with entries `Naive`, `Oracle`, `SIMEX`.

## Examples

``` r
method_colors()
#>     Naive    Oracle     SIMEX 
#> "#e07a5f" "#3d405b" "#81b29a" 
```
