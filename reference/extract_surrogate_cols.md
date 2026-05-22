# Identify surrogate-measurement columns in a data frame

Returns the names of columns matching a regex pattern. The default
`"^W[0-9]+$"` matches `W1`, `W2`, ... and corresponds to the convention
used by
[`generate_aft_data()`](https://qihuangzhang.github.io/aftsplex/reference/generate_aft_data.md)
and the manuscript notation.

## Usage

``` r
extract_surrogate_cols(data, pattern = "^W[0-9]+$")
```

## Arguments

- data:

  A data frame.

- pattern:

  A regular expression. Defaults to `"^W[0-9]+$"`.

## Value

A character vector of column names, in the order they appear in `data`.
Throws an error if no columns match.

## Examples

``` r
sim <- generate_aft_data(n = 50, n_val = 20, seed = 1)
extract_surrogate_cols(sim$survival)
#> [1] "W1" "W2" "W3"
```
