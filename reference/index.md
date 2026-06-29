# Package index

## Data-generating process

Helpers for the simulation studies used in the manuscript.

- [`f_true()`](https://qihuangzhang.github.io/aftsplex/reference/f_true.md)
  : True dose-response used by the simulation DGP
- [`build_sigma_u()`](https://qihuangzhang.github.io/aftsplex/reference/build_sigma_u.md)
  : Build a measurement-error covariance from marginal variances and a
  correlation
- [`sigma_u_for_reliability()`](https://qihuangzhang.github.io/aftsplex/reference/sigma_u_for_reliability.md)
  : Build a measurement-error covariance targeting a given reliability
- [`generate_aft_data()`](https://qihuangzhang.github.io/aftsplex/reference/generate_aft_data.md)
  : Generate one Monte Carlo replicate of the AFT-spline-SIMEX scenario

## Phase 1 — Measurement-error calibration

Multivariate ME model and the GLS combiner.

- [`fit_me_calibration()`](https://qihuangzhang.github.io/aftsplex/reference/fit_me_calibration.md)
  : Fit the Phase-1 multivariate ME calibration model
- [`gls_combine()`](https://qihuangzhang.github.io/aftsplex/reference/gls_combine.md)
  : GLS-combine surrogate columns into a single calibrated exposure
- [`extract_surrogate_cols()`](https://qihuangzhang.github.io/aftsplex/reference/extract_surrogate_cols.md)
  : Identify surrogate-measurement columns in a data frame

## Phase 2 — AFT-spline fitting and SIMEX

The estimator pipeline.

- [`fit_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/fit_aft_spline.md)
  : Fit an AFT model with a natural-cubic-spline exposure
- [`predict_curve()`](https://qihuangzhang.github.io/aftsplex/reference/predict_curve.md)
  : Predict the centred dose-response curve at a reference covariate
  setting
- [`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md)
  : SIMEX-correct the AFT-spline dose-response curve
- [`summary(`*`<simex_aft_spline>`*`)`](https://qihuangzhang.github.io/aftsplex/reference/summary.simex_aft_spline.md)
  : Summarise a SIMEX-corrected AFT-spline fit
- [`plot(`*`<simex_aft_spline>`*`)`](https://qihuangzhang.github.io/aftsplex/reference/plot.simex_aft_spline.md)
  : Plot a SIMEX-corrected AFT-spline fit

## Inference — two-stage bootstrap

- [`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)
  : Two-stage nonparametric bootstrap for the SIMEX-corrected curve
- [`summary(`*`<two_stage_bootstrap>`*`)`](https://qihuangzhang.github.io/aftsplex/reference/summary.two_stage_bootstrap.md)
  : Summarise a two-stage bootstrap SIMEX dose-response fit

## Evaluation and visualisation

- [`compute_ise()`](https://qihuangzhang.github.io/aftsplex/reference/compute_ise.md)
  : Integrated squared error of a fitted curve against the truth
- [`plot_curves()`](https://qihuangzhang.github.io/aftsplex/reference/plot_curves.md)
  : Overlay fitted dose-response curves with the truth and an optional
  CI
- [`stagill_palette()`](https://qihuangzhang.github.io/aftsplex/reference/stagill_palette.md)
  : StaGILL lab palette for this project
- [`method_colors()`](https://qihuangzhang.github.io/aftsplex/reference/method_colors.md)
  : Semantic method-to-colour mapping for Naive / Oracle / SIMEX
  comparisons
