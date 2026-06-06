# Coverage diagnostics: bootstrap SE versus empirical SE

The two-stage bootstrap pointwise 95% interval for the SIMEX-corrected
curve attains empirical coverage around 85% on the simulated
diminishing- returns curve — short of the nominal 95%. This vignette
reproduces the diagnostic from the manuscript appendix that pinpoints
the cause: the bootstrap recovers the sampling variance honestly, and
the coverage gap is bias-centering, not variance-underestimation.

``` r

library(aftsplex)

N_REP    <- 10L     # outer Monte Carlo replicates (manuscript uses 100)
R_BOOT   <- 30L     # inner bootstrap replicates (manuscript uses 200)
B_INNER  <- 10L
N_MAIN   <- 1500
N_VAL    <- 400
DF       <- 4
V_REF    <- c(V1 = 30, V2 = 30, V3 = 0, V4 = 0)

x_grid <- seq(qnorm(0.02, 10, sqrt(5)),
              qnorm(0.98, 10, sqrt(5)),
              length.out = 60)
truth  <- f_true(x_grid) - f_true(x_grid[1])
```

## One replicate

``` r

run_one <- function(seed) {
  sim <- generate_aft_data(n = N_MAIN, n_val = N_VAL, seed = seed)
  boot <- two_stage_bootstrap(
    survival   = sim$survival,
    validation = sim$validation,
    x_var      = "X_true",
    covariates = names(V_REF), v_ref = V_REF,
    df         = DF,
    lambda     = c(0.5, 1, 1.5, 2),
    B = B_INNER, R = R_BOOT,
    x_grid = x_grid
  )
  list(
    f_hat   = boot$f_hat,
    width   = boot$upper - boot$lower,
    covered = as.integer(truth >= boot$lower & truth <= boot$upper)
  )
}
```

## Monte Carlo

``` r

t0 <- Sys.time()
res <- lapply(seq_len(N_REP), function(i) run_one(20260601L + i))
elapsed <- Sys.time() - t0
```

``` r

elapsed
#> Time difference of 3.118939 mins
```

## Bootstrap SE versus empirical SE

For each grid point, the bootstrap CI width implies a
normal-approximation SE of `width / (2 * 1.96)`. The empirical SE is the
Monte Carlo SD of the single-fit `f_hat` across replicates.

``` r

width_mat <- do.call(cbind, lapply(res, `[[`, "width"))
fhat_mat  <- do.call(cbind, lapply(res, `[[`, "f_hat"))

se_boot <- rowMeans(width_mat) / (2 * qnorm(0.975))
se_emp  <- apply(fhat_mat, 1, sd)

interior <- x_grid >= quantile(x_grid, 0.05) &
            x_grid <= quantile(x_grid, 0.95)

data.frame(
  region  = c("Full grid", "Interior 5-95"),
  SE_boot = round(c(mean(se_boot),           mean(se_boot[interior])), 3),
  SE_emp  = round(c(mean(se_emp),            mean(se_emp[interior])), 3),
  ratio   = round(c(mean(se_boot) / mean(se_emp),
                    mean(se_boot[interior]) / mean(se_emp[interior])), 3)
)
#>          region SE_boot SE_emp ratio
#> 1     Full grid   0.098  0.090 1.079
#> 2 Interior 5-95   0.101  0.095 1.065
```

In the full coverage study (`N_REP = 100`, `R_BOOT = 200`) the ratio is
`0.999` (full grid) and `1.000` (interior) — the bootstrap is honest.

## Bias as the residual explanation

``` r

bias <- rowMeans(fhat_mat) - truth
data.frame(
  region   = c("Full grid", "Interior 5-95"),
  abs_bias = round(c(mean(abs(bias)),           mean(abs(bias[interior]))), 3),
  abs_SE   = round(c(mean(se_emp),              mean(se_emp[interior])), 3),
  ratio    = round(c(mean(abs(bias)) / mean(se_emp),
                     mean(abs(bias[interior])) / mean(se_emp[interior])), 3)
)
#>          region abs_bias abs_SE ratio
#> 1     Full grid    0.039  0.090 0.436
#> 2 Interior 5-95    0.041  0.095 0.434
```

When the bias-to-SE ratio sits well above zero, the coverage gap is a
centering problem. The manuscript decomposes this further into a
df-spline-smoothing portion (shared with the Oracle estimator) and a
SIMEX-specific portion (~0.02), and shows that higher spline `df`
reduces the smoothing portion modestly but does not eliminate it.
