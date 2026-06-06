# Fit the Phase-1 multivariate ME calibration model

Fits per-equation OLS of each surrogate column on the validation truth:
`W_j = alpha_0j + alpha_1j * X + e_j`, with `e ~ N_J(0, Sigma_e)`. Under
shared regressors this is equivalent to seemingly-unrelated regression
(SUR) and admits closed-form ML.

## Usage

``` r
fit_me_calibration(
  validation,
  x_var = "X_true",
  surrogate_pattern = "^W[0-9]+$"
)
```

## Arguments

- validation:

  Validation data frame containing the truth column named in `x_var` and
  the surrogate columns matched by `surrogate_pattern`.

- x_var:

  Name of the truth column. Default `"X_true"`.

- surrogate_pattern:

  Regex used to identify surrogate columns; see
  [`extract_surrogate_cols()`](https://qihuangzhang.github.io/aftsplex/reference/extract_surrogate_cols.md).
  Default `"^W[0-9]+$"`.

## Value

A list with elements

- `alpha0` numeric vector of intercepts (length `J`),

- `alpha1` numeric vector of slopes (length `J`),

- `Sigma_e` `J x J` residual covariance,

- `W_cols` character vector of surrogate column names used,

- `n_val` number of validation rows.

## Examples

``` r
sim <- generate_aft_data(n = 100, n_val = 200, seed = 1)
fit_me_calibration(sim$validation)
#> $alpha0
#> [1] 0.2168513 0.8582034 0.3415081
#> 
#> $alpha1
#> [1] 0.9775446 0.9069137 0.9746490
#> 
#> $Sigma_e
#>           [,1]      [,2]      [,3]
#> [1,] 2.1941712 0.8586757 0.9393987
#> [2,] 0.8586757 2.4148451 0.9587875
#> [3,] 0.9393987 0.9587875 2.7902613
#> 
#> $W_cols
#> [1] "W1" "W2" "W3"
#> 
#> $n_val
#> [1] 200
#> 
```
