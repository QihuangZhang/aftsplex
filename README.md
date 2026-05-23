# aftsplex

<!-- badges: start -->
[![R-CMD-check](https://github.com/QihuangZhang/aftsplex/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/QihuangZhang/aftsplex/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/QihuangZhang/aftsplex/actions/workflows/pkgdown.yaml/badge.svg)](https://qihuangzhang.github.io/aftsplex/)
<!-- badges: end -->

Documentation: <https://qihuangzhang.github.io/aftsplex/>

Accelerated-failure-time (AFT) models with a natural-cubic-spline
dose-response, corrected for additive measurement error in a multivariate
surrogate via a generalised-least-squares (GLS) combiner and SIMEX, with a
two-stage nonparametric bootstrap that propagates uncertainty from both the
validation (Phase-1) and main-study (Phase-2) samples.

## Installation

```r
# install.packages("remotes")
remotes::install_github("QihuangZhang/aftsplex")
```

## Quick start

```r
library(aftsplex)

sim <- generate_aft_data(n = 2000, n_val = 500, seed = 1)

boot <- two_stage_bootstrap(
  survival   = sim$survival,
  validation = sim$validation,
  x_var      = "X_true",
  covariates = c("V1", "V2", "V3", "V4"),
  v_ref      = c(V1 = 30, V2 = 30, V3 = 0, V4 = 0),
  df         = 4,
  lambda     = seq(0.5, 2, 0.5),
  B          = 20,
  R          = 50
)

plot(boot$x_grid, boot$f_hat, type = "l",
     ylim = range(boot$lower, boot$upper),
     xlab = "X", ylab = "Centred linear predictor")
lines(boot$x_grid, boot$lower, lty = 2)
lines(boot$x_grid, boot$upper, lty = 2)
```

See `vignette("quickstart", package = "aftsplex")` for a fuller walk-through,
`vignette("simulation-study")` for the Monte Carlo bias and ISE comparison
across estimators, and `vignette("sensitivity-df")` for guidance on choosing
the spline degrees of freedom.

## Method

`aftsplex` implements a two-phase estimator:

1. **Phase 1.** Multivariate linear measurement-error calibration on an
   external validation sample with paired (truth, surrogate) observations:
   `log W*_ij = alpha_0j + alpha_1j * log W_i + e_ij`. A GLS combiner
   weights the back-transformed surrogates into a single calibrated
   exposure with conditional variance.
2. **Phase 2.** Natural-cubic-spline AFT model on the calibrated exposure,
   corrected for residual attenuation via SIMEX (Cook and Stefanski, 1994).

Inference uses a nested two-stage nonparametric bootstrap that resamples
the validation and main-study samples jointly, addressing the
plug-in-variance critique that applies whenever the validation sample is
small relative to the main study.

## References

- Carroll, R. J., Kuchenhoff, H., Lombard, F., & Stefanski, L. A. (1996).
  Asymptotics for the SIMEX estimator in nonlinear measurement-error models.
  *Journal of the American Statistical Association*, 91, 242-250.
- Cook, J. R., & Stefanski, L. A. (1994). Simulation-extrapolation
  estimation in parametric measurement-error models.
  *Journal of the American Statistical Association*, 89, 1314-1328.

## License

MIT (c) 2026 Qihuang Zhang.
