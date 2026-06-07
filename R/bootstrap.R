#' Two-stage nonparametric bootstrap for the SIMEX-corrected curve
#'
#' For each replicate `r`:
#'   1. resample validation participants;
#'   2. refit the multivariate ME calibration via [fit_me_calibration()];
#'   3. resample main-study participants;
#'   4. form the bootstrap-combined surrogate via [gls_combine()] under the
#'      refit calibration;
#'   5. run [simex_aft_spline()] with the refit `sigma_w_sq` and store the
#'      corrected curve.
#'
#' Pointwise quantiles across `r` give CIs that propagate both Phase-1 and
#' Phase-2 uncertainty - this addresses the plug-in-variance critique that
#' the conditional `sigma_w_sq` alone yields anti-conservative intervals
#' when the validation sample is small relative to the main study.
#'
#' @param survival Main-study data frame; must contain the columns named in
#'   `x_var` (if a single calibrated exposure already exists), or the
#'   surrogate columns matched by `surrogate_pattern`, plus the outcome,
#'   status, and covariate columns.
#' @param validation Validation data frame containing the truth column
#'   named in `x_var` and the surrogate columns matched by
#'   `surrogate_pattern`.
#' @param x_var Name of the exposure column. For simulation use the truth
#'   column name (`"X_true"`); for the standard workflow this is the same
#'   name used in both `validation` and the bootstrap-recombined
#'   `survival` data.
#' @param covariates Character vector of confounder column names.
#' @param v_ref Named numeric vector of reference covariate values matching
#'   `covariates`. `NULL` if no covariates.
#' @param surrogate_pattern Regex used to identify surrogate columns;
#'   default `"^W[0-9]+$"`.
#' @param df,knots,lambda,B,dist Passed to [simex_aft_spline()].
#' @param outcome_var,status_var Passed to [simex_aft_spline()].
#' @param R Number of outer bootstrap replicates.
#' @param x_grid Numeric exposure grid for the curve. If `NULL`, taken as
#'   100 equally spaced points between the 2nd and 98th percentile of the
#'   full-sample combined surrogate.
#' @param nested Logical; if `TRUE` (default) run the full two-stage
#'   bootstrap that resamples the validation sample and refits the Phase-1
#'   calibration on every replicate. If `FALSE`, hold the calibration fixed
#'   at the full-sample fit and resample only the main study (single-stage),
#'   giving the conditional-variance interval. Running both and comparing
#'   `mean(upper - lower)` quantifies the Phase-1 contribution to the width.
#' @param x_ref Optional scalar reference exposure at which to anchor the
#'   curves (and every bootstrap replicate) before taking quantiles, so the
#'   interval is read relative to `x_ref`. `NULL` (default) anchors at
#'   `x_grid[1]`.
#' @param support_probs Optional length-2 probabilities defining the
#'   trustworthy exposure support as quantiles of the full-sample combined
#'   surrogate. When set, an `in_support` vector is returned and a one-shot
#'   `warning()` flags out-of-support grid points. `NULL` (default) disables.
#' @param workers Integer number of parallel workers for the outer bootstrap
#'   loop. `1L` (default) runs serially and is the canonical reproducible path:
#'   it draws from R's global RNG stream, so a given `set.seed()` gives
#'   bit-for-bit identical results. When `> 1`, replicates run via
#'   `future.apply::future_lapply()` on a `multisession` plan with per-task
#'   L'Ecuyer streams (`future.seed = TRUE`); these are reproducible and
#'   identical *across worker counts* (e.g. `workers = 2` matches `workers = 3`),
#'   but use independent streams and therefore **differ from the serial
#'   (`workers = 1`) result for the same seed** by design. Requires the `future`
#'   and `future.apply` packages.
#' @param verbose Print progress bar? (Serial path only.)
#'
#' @return A list with elements
#'   * `x_grid` exposure grid,
#'   * `f_hat` single-fit SIMEX point estimate (full-sample, no resampling),
#'   * `bs_mean` mean of the bootstrap curves,
#'   * `curves` matrix of bootstrap curves (`length(x_grid)` rows by `R` cols),
#'   * `median, lower, upper` percentile-CI summaries (2.5%, 50%, 97.5%),
#'   * `lower_basic, upper_basic` reverse-percentile (basic) CI bounds,
#'   * `f_bc` bias-corrected point estimate `2*f_hat - bs_mean`,
#'   * `sigma_w_sq_boot, omega_boot` per-replicate Phase-1 summaries,
#'   * `in_support` logical support flag (or `NULL`),
#'   * `R_effective` number of replicates that produced a curve.
#'
#' The percentile CI is the primary interval; the basic (reverse-percentile)
#' CI is kept for completeness but does not improve coverage on simulated
#' diminishing-returns curves because the dominant residual bias is from
#' spline smoothing rather than SIMEX.
#'
#' @examples
#' \donttest{
#' sim <- generate_aft_data(n = 500, n_val = 200, seed = 1)
#' boot <- two_stage_bootstrap(
#'   survival   = sim$survival,
#'   validation = sim$validation,
#'   x_var      = "X_true",
#'   covariates = c("V1","V2","V3","V4"),
#'   v_ref      = c(V1=30, V2=30, V3=0, V4=0),
#'   lambda     = c(0.5, 1, 1.5, 2), B = 5, R = 10
#' )
#' head(boot$f_hat)
#'
#' # Phase-1 contribution to CI width: compare double vs single-stage.
#' single <- two_stage_bootstrap(
#'   survival = sim$survival, validation = sim$validation, x_var = "X_true",
#'   covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30, V2=30, V3=0, V4=0),
#'   lambda = c(0.5, 1, 1.5, 2), B = 5, R = 10, nested = FALSE
#' )
#' mean(boot$upper - boot$lower) / mean(single$upper - single$lower)
#' }
#'
#' @export
two_stage_bootstrap <- function(survival, validation, x_var = "X_true",
                                covariates = character(0), v_ref = NULL,
                                surrogate_pattern = "^W[0-9]+$",
                                df = 4, knots = NULL,
                                outcome_var = "T_obs", status_var = "delta",
                                lambda = c(0.5, 1, 1.5, 2), B = 50,
                                R = 200, x_grid = NULL,
                                nested = TRUE, x_ref = NULL,
                                support_probs = NULL, workers = 1L,
                                dist = "lognormal", verbose = FALSE) {
  if (length(workers) != 1L || !is.finite(workers) || workers < 1L) {
    stop("'workers' must be a positive integer (1 = serial).", call. = FALSE)
  }
  if (!is.null(v_ref) && length(v_ref) != length(covariates)) {
    stop("'v_ref' has length ", length(v_ref), " but 'covariates' has length ",
         length(covariates), "; they must match.", call. = FALSE)
  }

  W_cols <- extract_surrogate_cols(survival, surrogate_pattern)
  W_mat_full <- as.matrix(survival[, W_cols, drop = FALSE])

  combo_full <- tryCatch({
    cal_full <- fit_me_calibration(validation, x_var = x_var,
                                   surrogate_pattern = surrogate_pattern)
    gls_combine(W_mat_full, cal_full)
  }, error = function(e) {
    stop("Full-sample Phase-1 calibration failed: ", conditionMessage(e),
         ". This is usually a singular calibration from collinear or ",
         "constant surrogate columns; check the '", surrogate_pattern,
         "' columns in 'validation'.", call. = FALSE)
  })
  W_bar_full <- combo_full$W_bar

  if (is.null(x_grid)) {
    x_grid <- seq(quantile(W_bar_full, 0.02),
                  quantile(W_bar_full, 0.98), length.out = 100)
  }

  data_full <- survival
  data_full$W_bar <- W_bar_full
  sim_full <- tryCatch(
    suppressWarnings(simex_aft_spline(
      data = data_full, x_var = "W_bar",
      sigma_w_sq = combo_full$sigma_w_sq,
      covariates = covariates, v_ref = v_ref,
      df = df, knots = knots,
      outcome_var = outcome_var, status_var = status_var,
      lambda = lambda, B = B,
      dist = dist, x_grid = x_grid
    )),
    error = function(e) NULL
  )
  f_hat <- if (!is.null(sim_full)) sim_full$curve_simex else rep(NA_real_, length(x_grid))

  curves <- matrix(NA_real_, length(x_grid), R)
  sigma_w_sq_boot <- rep(NA_real_, R)
  omega_boot <- matrix(NA_real_, length(W_cols), R)

  # One bootstrap replicate. Self-contained (own resamples, own Phase-1 fit
  # when `nested`), so it is the parallel unit. Calls below keep the same RNG
  # call order as the original loop, so the serial path is bit-for-bit unchanged.
  run_rep <- function(r) {
    # The whole replicate - resampling, Phase-1 refit, GLS combine and SIMEX -
    # is wrapped so that a degenerate resample (e.g. a singular calibration or
    # a non-convergent fit) drops just that replicate to NA instead of
    # aborting the entire bootstrap. Per-replicate warnings are muted; the
    # caller reports a single summary of how many replicates were lost.
    tryCatch(
      suppressWarnings({
        if (nested) {
          v_idx <- sample.int(nrow(validation), replace = TRUE)
          cal_r <- fit_me_calibration(validation[v_idx, , drop = FALSE],
                                      x_var = x_var,
                                      surrogate_pattern = surrogate_pattern)
        } else {
          cal_r <- cal_full
        }

        s_idx <- sample.int(nrow(survival), replace = TRUE)
        data_r <- survival[s_idx, , drop = FALSE]

        W_mat_r <- as.matrix(data_r[, W_cols, drop = FALSE])
        combo_r <- gls_combine(W_mat_r, cal_r)
        data_r$W_bar <- combo_r$W_bar

        sim_r <- simex_aft_spline(
          data = data_r, x_var = "W_bar",
          sigma_w_sq = combo_r$sigma_w_sq,
          covariates = covariates, v_ref = v_ref,
          df = df, knots = knots,
          outcome_var = outcome_var, status_var = status_var,
          lambda = lambda, B = B,
          dist = dist, x_grid = x_grid
        )
        list(curve = sim_r$curve_simex,
             sigma_w_sq = combo_r$sigma_w_sq,
             omega = combo_r$omega)
      }),
      error = function(e) NULL
    )
  }

  if (workers > 1L) {
    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop("workers > 1 requires the 'future' and 'future.apply' packages.")
    }
    oplan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(oplan), add = TRUE)
    res_list <- future.apply::future_lapply(seq_len(R), run_rep,
                                            future.seed = TRUE)
  } else {
    pb <- if (verbose) txtProgressBar(min = 0, max = R, style = 3) else NULL
    res_list <- vector("list", R)
    for (r in seq_len(R)) {
      # Single-bracket assignment: a failed replicate returns NULL, and
      # `res_list[[r]] <- NULL` would *delete* the element (shrinking the list
      # below R and breaking the `ran`/`R_effective` bookkeeping below). Wrapping
      # in list() stores an explicit NULL so the length stays R.
      res_list[r] <- list(run_rep(r))
      if (!is.null(pb)) setTxtProgressBar(pb, r)
    }
    if (!is.null(pb)) { close(pb); message("") }
  }

  ran <- !vapply(res_list, is.null, logical(1))
  for (r in seq_len(R)) {
    if (ran[r]) {
      rr <- res_list[[r]]
      curves[, r]        <- rr$curve
      sigma_w_sq_boot[r] <- rr$sigma_w_sq
      omega_boot[, r]    <- rr$omega
    }
  }
  # A replicate is "effective" only if it yielded at least one finite curve
  # value; a replicate can run (valid calibration) yet still produce an all-NA
  # curve when every perturbed SIMEX fit fails to converge.
  R_effective <- sum(colSums(is.finite(curves)) > 0L)

  # Tell the user what went wrong rather than returning silent NAs (the curves
  # of failed replicates are all NA and quietly dropped by na.rm downstream).
  n_failed <- R - R_effective
  if (R_effective == 0L) {
    warning(sprintf(
      paste0("All %d bootstrap replicates failed to produce a curve; every ",
             "confidence-interval summary is NA. This usually indicates too ",
             "few events, heavy censoring, or a singular calibration. Inspect ",
             "the data, reduce 'df', or increase the sample size."), R),
      call. = FALSE)
  } else if (n_failed > 0L) {
    warning(sprintf(
      paste0("%d of %d bootstrap replicates failed and were dropped (NA); ",
             "intervals are based on the remaining %d. Treat the interval ",
             "with caution if this fraction is large."),
      n_failed, R, R_effective), call. = FALSE)
  }
  grid_na <- rowSums(!is.na(curves)) == 0L
  if (R_effective > 0L && any(grid_na)) {
    warning(sprintf(
      paste0("%d of %d grid point(s) have no finite bootstrap value (CI is NA ",
             "there); these are typically spline extrapolations beyond the ",
             "exposure support."),
      sum(grid_na), length(x_grid)), call. = FALSE)
  }
  if (all(is.na(f_hat))) {
    warning("The full-sample SIMEX point estimate ('f_hat') is entirely NA ",
            "(the full-sample fit failed or could not be extrapolated); ",
            "only bootstrap summaries are available.", call. = FALSE)
  }

  if (!is.null(x_ref)) {
    f_hat <- recenter_curve(f_hat, x_grid, x_ref)
    for (r in seq_len(R)) {
      if (!all(is.na(curves[, r]))) {
        curves[, r] <- recenter_curve(curves[, r], x_grid, x_ref)
      }
    }
  }

  in_support <- support_flag(x_grid, W_bar_full, support_probs)

  # Robust pointwise summaries: a grid row with no finite values yields NA
  # instead of erroring (quantile() on an empty vector) or NaN (rowMeans).
  safe_quantile <- function(v, p) {
    v <- v[is.finite(v)]
    if (length(v) == 0L) NA_real_ else unname(quantile(v, probs = p))
  }
  safe_median <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0L) NA_real_ else median(v)
  }
  q025 <- apply(curves, 1, safe_quantile, p = 0.025)
  q975 <- apply(curves, 1, safe_quantile, p = 0.975)
  bs_mean <- rowMeans(curves, na.rm = TRUE)
  bs_mean[is.nan(bs_mean)] <- NA_real_

  list(
    x_grid      = x_grid,
    f_hat       = f_hat,
    bs_mean     = bs_mean,
    curves      = curves,
    median      = apply(curves, 1, safe_median),
    lower       = q025,
    upper       = q975,
    lower_basic = 2 * f_hat - q975,
    upper_basic = 2 * f_hat - q025,
    f_bc        = 2 * f_hat - bs_mean,
    sigma_w_sq_boot = sigma_w_sq_boot,
    omega_boot      = omega_boot,
    in_support      = in_support,
    R_effective     = R_effective
  )
}
