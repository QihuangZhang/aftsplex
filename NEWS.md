# aftsplex 0.3.1

Robustness fixes surfaced by an adversarial stress-test of 0.3.0. All changes are
input guards and a bug fix; valid-input behaviour and the byte-identical
`generate_aft_data()` baseline are unchanged.

* `two_stage_bootstrap()` no longer crashes on the serial path (`workers = 1`,
  the default) when bootstrap replicates fail. A failed replicate returned `NULL`
  and `list[[r]] <- NULL` deleted the element, shrinking the result list and
  triggering "missing value where TRUE/FALSE needed"; the documented
  drop-and-warn behaviour now works on both the serial and parallel paths.
* `fit_me_calibration()`, `gls_combine()`, `compute_ise()`,
  `sigma_u_for_reliability()`, `build_sigma_u()`, and `generate_aft_data()` now
  validate their inputs and raise informative errors instead of silently
  propagating `NaN`/`Inf` or returning a meaningless result: constant/too-small
  validation samples, non-finite or uninformative (zero-slope) surrogates,
  singular surrogate covariances, unsorted or length-mismatched `compute_ise()`
  arguments, non-positive `var_x`/`ratios`/`sigma_u_sq`, out-of-range `rho`, and
  negative `x_sd`.
* Clarified that the serial (`workers = 1`) bootstrap is the canonical
  bit-for-bit reproducible path; parallel runs use independent L'Ecuyer streams
  that are reproducible across worker counts but differ from the serial result
  for the same seed by design.

# aftsplex 0.3.0

Realistic-scale data-generating process, two new simulation diagnostics, and a
robustness pass on the bootstrap/SIMEX failure modes.

* `generate_aft_data()` gains `x_mean`, `x_sd`, and `f_args` so the latent
  exposure and dose-response can be expressed on a meaningful scale (e.g. daily
  light-physical-activity minutes, mean ~300 with a 150-minute threshold). The
  defaults reproduce the previous arbitrary-unit output bit-for-bit.
* New exported `sigma_u_for_reliability()` returns a measurement-error
  covariance targeting a given mean per-surrogate reliability, for sweeping
  measurement quality (e.g. 0.5 / 0.65 / 0.8).
* The `simulation-study` vignette gains a "daily-minutes scale" section with a
  reliability sweep (Oracle / Naive / SIMEX ISE across reliability) and a
  per-surrogate ladder (each surrogate alone vs the GLS combiner via the
  existing `surrogate_pattern` argument; no engine change).
* Robustness pass (a user hit an all-NA bootstrap curve with no message, then
  an opaque error at plotting). `two_stage_bootstrap()` now isolates each
  replicate so a degenerate resample drops to NA instead of aborting the run,
  and emits a single summary warning (replicates dropped, grid points with no
  finite value, all-NA point estimate, or all replicates failed). The
  full-sample calibration gives a clear message on collinear/constant
  surrogates instead of a raw Lapack error, and pointwise summaries are NA-safe.
* `simex_aft_spline()` returns NA at grid points with too few converged
  perturbed fits (one summary warning) instead of erroring with
  "0 (non-NA) cases"; a failed naive fit now gives an informative error.
* `plot_curves()` drops entirely-NA curves with a warning and stops with a
  clear message when nothing is left to plot, instead of an empty panel.
* New "Troubleshooting" article documenting every warning/error message, and a
  "coverage diagnostics" article.

# aftsplex 0.2.0

Revision round addressing reviewer feedback. All additions are backward-compatible
via defaults; existing output and the SIMEX regression are unchanged.

* `two_stage_bootstrap()` gains `nested` to toggle between the full two-stage
  (double) bootstrap and a single-stage interval that holds the Phase-1
  calibration fixed, so the Phase-1 contribution to CI width can be quantified.
* `two_stage_bootstrap()` and `simex_aft_spline()` gain `x_ref` to anchor the
  centred curves (and every bootstrap replicate) at a clinically meaningful
  reference exposure before quantiles are taken.
* `two_stage_bootstrap()` and `simex_aft_spline()` gain `support_probs`, which
  flags grid points outside the exposure support (returned as `in_support`) and
  warns that they are spline extrapolations rather than fits.
* `plot_curves()` and `plot.simex_aft_spline()` gain `time_ratio` to plot the
  exponentiated (time-ratio) curve with a reference line at 1, and render
  out-of-support segments dashed with a lighter ribbon.
* `two_stage_bootstrap()` gains `workers` for parallel execution of the outer
  bootstrap loop via `future.apply` (per-task L'Ecuyer streams; serial default
  unchanged). Adds `future` and `future.apply` to Suggests.

# aftsplex 0.1.0

* Initial release.
