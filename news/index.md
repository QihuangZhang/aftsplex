# Changelog

## aftsplex 0.2.0

Revision round addressing reviewer feedback. All additions are
backward-compatible via defaults; existing output and the SIMEX
regression are unchanged.

- [`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)
  gains `nested` to toggle between the full two-stage (double) bootstrap
  and a single-stage interval that holds the Phase-1 calibration fixed,
  so the Phase-1 contribution to CI width can be quantified.
- [`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)
  and
  [`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md)
  gain `x_ref` to anchor the centred curves (and every bootstrap
  replicate) at a clinically meaningful reference exposure before
  quantiles are taken.
- [`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)
  and
  [`simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/simex_aft_spline.md)
  gain `support_probs`, which flags grid points outside the exposure
  support (returned as `in_support`) and warns that they are spline
  extrapolations rather than fits.
- [`plot_curves()`](https://qihuangzhang.github.io/aftsplex/reference/plot_curves.md)
  and
  [`plot.simex_aft_spline()`](https://qihuangzhang.github.io/aftsplex/reference/plot.simex_aft_spline.md)
  gain `time_ratio` to plot the exponentiated (time-ratio) curve with a
  reference line at 1, and render out-of-support segments dashed with a
  lighter ribbon.
- [`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)
  gains `workers` for parallel execution of the outer bootstrap loop via
  `future.apply` (per-task L’Ecuyer streams; serial default unchanged).
  Adds `future` and `future.apply` to Suggests.

## aftsplex 0.1.0

- Initial release.
