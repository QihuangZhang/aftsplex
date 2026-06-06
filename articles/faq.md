# Troubleshooting: warnings and error messages

`aftsplex` is designed to **fail loudly and specifically**: when a fit
cannot be produced it tells you what happened and where, rather than
returning a silent `NA` that surfaces later as a confusing plotting
error. This page lists the messages you may see and what to check next.

Most issues trace back to a small set of causes:

- **too few events** (heavy censoring) for the spline degrees of
  freedom,
- a **singular calibration** from collinear or constant surrogate
  columns,
- an **`x_grid` that runs past the exposure support**, where the spline
  only extrapolates,
- a missing optional package.

## Warnings

Warnings mean the run completed but part of the result is `NA` or should
be read with care. The object is still returned, so you can inspect it.

| Message (pattern) | What it means | What to check / do |
|----|----|----|
| `N of R bootstrap replicates failed and were dropped …` | Some replicates hit a degenerate resample (non-convergent fit or singular calibration) and were skipped. The interval uses the survivors. | Fine if `N` is a small fraction. If large, increase the sample size or number of events, or lower `df`. Inspect `boot$R_effective`. |
| `All R bootstrap replicates failed to produce a curve …` | No replicate yielded a usable curve, so every CI summary is `NA`. | Almost always too few events, heavy censoring, or a singular calibration. Check the event count (`sum(delta)`), reduce `df`, increase `n`, and verify the surrogate columns vary. |
| `N of M grid point(s) have no finite bootstrap value …` | The CI is `NA` at some exposures even though other replicates succeeded — typically the grid extends past the spline’s support. | Trim `x_grid` to the bulk of the exposure, or set `support_probs` to flag (rather than hide) the extrapolated region. |
| `The full-sample SIMEX point estimate ('f_hat') is entirely NA …` | The point-estimate fit on the full sample failed or could not be extrapolated; only the bootstrap summaries are available. | Same causes as above — check events and `df`. The bootstrap median/quantiles may still be usable. |
| `SIMEX extrapolation undefined at N of M grid point(s) …` | At those grid points fewer than three of the perturbed-`lambda` fits converged, so the quadratic extrapolation to `lambda = -1` is not defined and returns `NA`. | Reduce `df`, lower the largest `lambda`, or increase `B`. Often co-occurs with heavy censoring. |
| `Curve(s) … are entirely NA and were dropped from the plot …` | [`plot_curves()`](https://qihuangzhang.github.io/aftsplex/reference/plot_curves.md) received a fitted curve that is all `NA` and omitted it. | Check `R_effective` and the warnings above for that fit before plotting. |
| `Confidence interval is entirely NA; omitting the ribbon` | The supplied `ci$lower`/`ci$upper` are all `NA`, so no ribbon is drawn. | The bootstrap produced no finite interval; see the all-replicates-failed guidance. |
| `N of M grid point(s) lie outside the support […]` | Informational: with `support_probs` set, these grid points are spline extrapolations beyond the trustworthy exposure range. | Expected when you deliberately show results past the support. They are drawn dashed; read them as extrapolations. |

## Errors

Errors stop execution because the result would be meaningless. They are
raised with a clear cause instead of a low-level message (e.g. a raw
`Lapack` singularity).

| Message (pattern) | What it means | What to check / do |
|----|----|----|
| `Full-sample Phase-1 calibration failed: … collinear or constant surrogate columns` | The calibration matrix is singular — two surrogate columns are (near-)identical, or the validation truth has no variation. | Drop duplicated surrogate columns, confirm each varies, and check the validation `x_var` is not constant. |
| `SIMEX naive AFT fit failed: …` | The underlying `survreg` fit could not be estimated at all (e.g. no events, or a rank-deficient spline basis). | Ensure there are enough events, reduce `df`, or simplify the model. |
| `Nothing to plot: every curve is entirely NA` | Every curve passed to [`plot_curves()`](https://qihuangzhang.github.io/aftsplex/reference/plot_curves.md) is `NA`, so there is nothing to draw. | Inspect the SIMEX/bootstrap output — `R_effective` and the failure warnings — before plotting. |
| `workers > 1 requires the 'future' and 'future.apply' packages` | Parallel execution was requested but the packages are not installed. | `install.packages(c("future", "future.apply"))`, or run serially with `workers = 1L` (the default). |
| `Package 'ggplot2' is required for plot_curves()` | The plotting helper needs `ggplot2`. | `install.packages("ggplot2")`. |
| `No surrogate columns matched pattern '…'` | No columns matched `surrogate_pattern` (default `^W[0-9]+$`). | Pass the correct `surrogate_pattern`, or rename the surrogate columns to `W1`, `W2`, … |

## Checklist: diagnosing a failed bootstrap

When
[`two_stage_bootstrap()`](https://qihuangzhang.github.io/aftsplex/reference/two_stage_bootstrap.md)
warns or returns `NA` summaries:

1.  **`boot$R_effective`** — how many replicates produced a usable
    curve. If this is much smaller than `R`, the fits are failing.
2.  **Event count** — `sum(survival$delta)`. A natural-cubic-spline AFT
    needs enough events to support `df`; with few events, lower `df`.
3.  **Surrogate columns** — confirm the `^W[0-9]+$` columns are not
    duplicated or constant (a singular calibration aborts with a clear
    error).
4.  **Grid vs. support** — if only the edges of the curve are `NA`, the
    grid runs past the support; trim `x_grid` or set `support_probs` to
    flag it.
5.  **`in_support`** — the returned flag marks which grid points are
    within the trustworthy range.

See the [Quick
start](https://qihuangzhang.github.io/aftsplex/articles/quickstart.md)
article for runnable examples of `support_probs`, `x_ref`, and the other
0.2.0 options.
