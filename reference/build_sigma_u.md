# Build a measurement-error covariance from marginal variances and a correlation

Assembles a `J x J` covariance with the supplied marginal variances on
the diagonal and a common pairwise correlation `rho` off the diagonal.
Used to configure the multivariate surrogate noise in
[`generate_aft_data()`](https://qihuangzhang.github.io/aftsplex/reference/generate_aft_data.md).

## Usage

``` r
build_sigma_u(sigma_u_sq = c(2, 2.5, 3), rho = 0.4)
```

## Arguments

- sigma_u_sq:

  Numeric vector of marginal variances (length `J`).

- rho:

  Common pairwise correlation in `[-1/(J-1), 1]`.

## Value

A `J x J` numeric matrix.

## Examples

``` r
build_sigma_u(c(2.0, 2.5, 3.0), rho = 0.4)
#>           [,1]      [,2]      [,3]
#> [1,] 2.0000000 0.8944272 0.9797959
#> [2,] 0.8944272 2.5000000 1.0954451
#> [3,] 0.9797959 1.0954451 3.0000000
```
