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
#' @param verbose Print progress bar?
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
#'   * `R_effective` number of replicates that produced a curve.
#'
#' Empirically (see `vignette("coverage-diagnostics")`) the percentile CI is
#' the primary interval; the basic CI is kept for completeness but does not
#' improve coverage on simulated diminishing-returns curves because the
#' dominant residual bias is from spline smoothing rather than SIMEX.
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
                                dist = "lognormal", verbose = FALSE) {
  W_cols <- extract_surrogate_cols(survival, surrogate_pattern)
  W_mat_full <- as.matrix(survival[, W_cols, drop = FALSE])

  cal_full <- fit_me_calibration(validation, x_var = x_var,
                                 surrogate_pattern = surrogate_pattern)
  combo_full <- gls_combine(W_mat_full, cal_full)
  W_bar_full <- combo_full$W_bar

  if (is.null(x_grid)) {
    x_grid <- seq(quantile(W_bar_full, 0.02),
                  quantile(W_bar_full, 0.98), length.out = 100)
  }

  data_full <- survival
  data_full$W_bar <- W_bar_full
  sim_full <- tryCatch(
    simex_aft_spline(
      data = data_full, x_var = "W_bar",
      sigma_w_sq = combo_full$sigma_w_sq,
      covariates = covariates, v_ref = v_ref,
      df = df, knots = knots,
      outcome_var = outcome_var, status_var = status_var,
      lambda = lambda, B = B,
      dist = dist, x_grid = x_grid
    ),
    error = function(e) NULL
  )
  f_hat <- if (!is.null(sim_full)) sim_full$curve_simex else rep(NA_real_, length(x_grid))

  curves <- matrix(NA_real_, length(x_grid), R)
  sigma_w_sq_boot <- numeric(R)
  omega_boot <- matrix(NA_real_, length(W_cols), R)

  pb <- if (verbose) txtProgressBar(min = 0, max = R, style = 3) else NULL

  for (r in seq_len(R)) {
    v_idx <- sample.int(nrow(validation), replace = TRUE)
    cal_r <- fit_me_calibration(validation[v_idx, , drop = FALSE],
                                x_var = x_var,
                                surrogate_pattern = surrogate_pattern)

    s_idx <- sample.int(nrow(survival), replace = TRUE)
    data_r <- survival[s_idx, , drop = FALSE]

    W_mat_r <- as.matrix(data_r[, W_cols, drop = FALSE])
    combo_r <- gls_combine(W_mat_r, cal_r)
    data_r$W_bar <- combo_r$W_bar

    sim_r <- tryCatch(
      simex_aft_spline(
        data = data_r, x_var = "W_bar",
        sigma_w_sq = combo_r$sigma_w_sq,
        covariates = covariates, v_ref = v_ref,
        df = df, knots = knots,
        outcome_var = outcome_var, status_var = status_var,
        lambda = lambda, B = B,
        dist = dist, x_grid = x_grid
      ),
      error = function(e) NULL
    )

    if (!is.null(sim_r)) {
      curves[, r] <- sim_r$curve_simex
      sigma_w_sq_boot[r] <- combo_r$sigma_w_sq
      omega_boot[, r] <- combo_r$omega
    }
    if (!is.null(pb)) setTxtProgressBar(pb, r)
  }
  if (!is.null(pb)) { close(pb); message("") }

  q025 <- apply(curves, 1, quantile, probs = 0.025, na.rm = TRUE)
  q975 <- apply(curves, 1, quantile, probs = 0.975, na.rm = TRUE)
  bs_mean <- rowMeans(curves, na.rm = TRUE)

  list(
    x_grid      = x_grid,
    f_hat       = f_hat,
    bs_mean     = bs_mean,
    curves      = curves,
    median      = apply(curves, 1, median, na.rm = TRUE),
    lower       = q025,
    upper       = q975,
    lower_basic = 2 * f_hat - q975,
    upper_basic = 2 * f_hat - q025,
    f_bc        = 2 * f_hat - bs_mean,
    sigma_w_sq_boot = sigma_w_sq_boot,
    omega_boot      = omega_boot,
    R_effective     = sum(!is.na(curves[1, ]))
  )
}
