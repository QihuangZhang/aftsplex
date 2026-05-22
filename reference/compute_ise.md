# Integrated squared error of a fitted curve against the truth

Trapezoidal-rule approximation of `integral of (f_hat - f_true)^2 dx`
over the exposure grid. The truth function is centred at the lower grid
boundary to match the convention used elsewhere in the package
(`truth(x_grid) - truth(x_grid[1])`).

## Usage

``` r
compute_ise(x_grid, truth_fn, fits)
```

## Arguments

- x_grid:

  Numeric vector of exposure values.

- truth_fn:

  A function such as
  [`f_true()`](https://qihuangzhang.github.io/aftsplex/reference/f_true.md)
  returning the (uncentred) true response at `x_grid`.

- fits:

  A named list of numeric vectors, each of length `length(x_grid)`,
  holding fitted centred curves.

## Value

A list with elements `truth` (the centred truth at `x_grid`) and `ise`
(named numeric vector of ISEs, one entry per element of `fits`).

## Examples

``` r
x_grid <- seq(6, 14, length.out = 50)
fits <- list(naive = numeric(length(x_grid)),
             oracle = f_true(x_grid) - f_true(x_grid[1]))
compute_ise(x_grid, f_true, fits)
#> $truth
#>  [1] 0.00000000 0.03793161 0.07346521 0.10675240 0.13793519 0.16714663
#>  [7] 0.19451133 0.22014606 0.24416018 0.26665615 0.28772993 0.30747144
#> [13] 0.32596490 0.34328922 0.35951831 0.37472141 0.38896338 0.40230497
#> [19] 0.41480312 0.42651115 0.43747900 0.44775347 0.45737840 0.46639485
#> [25] 0.47484128 0.48275373 0.49016596 0.49710960 0.50361426 0.50970771
#> [31] 0.51541593 0.52076328 0.52577257 0.53046518 0.53486113 0.53897916
#> [37] 0.54283686 0.54645068 0.54983603 0.55300737 0.55597821 0.55876124
#> [43] 0.56136833 0.56381060 0.56609847 0.56824170 0.57024944 0.57213025
#> [49] 0.57389216 0.57554268
#> 
#> $ise
#>   naive  oracle 
#> 1.60217 0.00000 
#> 
```
