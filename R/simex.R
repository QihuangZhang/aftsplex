#' SIMEX-correct the AFT-spline dose-response curve
#'
#' Implements simulation-extrapolation (Cook and Stefanski 1994; Carroll,
#' Kuchenhoff, Lombard and Stefanski 1996) for a natural-cubic-spline AFT
#' model fitted by [fit_aft_spline()]. For each `lambda` in the supplied
#' grid, draws `B` perturbed surrogates with additional noise
#' `sqrt(lambda * sigma_w_sq) * N(0, 1)`, refits, averages the centred
#' linear predictor across `B`, then extrapolates to `lambda = -1` via a
#' pointwise quadratic OLS.
#'
#' Inner `survreg` fits are wrapped in convergence handling: fits that
#' fail to converge within `maxiter` iterations are dropped from the
#' average via `na.rm = TRUE`, and the count is reported if `verbose = TRUE`.
#'
#' @param data A data frame containing the columns named in `x_var`,
#'   `covariates`, `outcome_var`, and `status_var`.
#' @param x_var Name of the exposure column.
#' @param sigma_w_sq Scalar conditional variance of the calibrated exposure
#'   (typically `gls_combine(...)$sigma_w_sq`).
#' @param covariates Character vector of confounder column names.
#' @param v_ref Named numeric vector of reference covariate values matching
#'   `covariates`. `NULL` if no covariates.
#' @param df,knots,dist Passed to [fit_aft_spline()].
#' @param outcome_var,status_var Passed to [fit_aft_spline()].
#' @param lambda Numeric grid of SIMEX lambda values (excluding 0; 0 is
#'   appended internally as the naive fit). Default `c(0.5, 1, 1.5, 2)`.
#' @param B Number of inner replicates per lambda.
#' @param x_grid Numeric grid of exposure values at which to evaluate the
#'   corrected curve. Defaults to 100 equally spaced points spanning
#'   `data[[x_var]]`.
#' @param verbose Print progress bar and non-convergence count?
#' @param maxiter `survreg` iteration cap (passed via
#'   [survival::survreg.control()]).
#'
#' @return A list with elements
#'   * `x_grid` the exposure grid,
#'   * `curve_simex` the centred SIMEX-corrected linear predictor,
#'   * `curves_lambda` matrix of centred lambda-fits (columns = lambda, including 0),
#'   * `lambda` the lambda values used (with 0 in column 1).
#'
#' @examples
#' \donttest{
#' sim <- generate_aft_data(n = 500, n_val = 200, seed = 1)
#' cal <- fit_me_calibration(sim$validation)
#' W   <- as.matrix(sim$survival[, cal$W_cols])
#' g   <- gls_combine(W, cal)
#' dat <- sim$survival
#' dat$W_bar <- g$W_bar
#' s <- simex_aft_spline(dat, x_var = "W_bar",
#'                       sigma_w_sq = g$sigma_w_sq,
#'                       covariates = c("V1","V2","V3","V4"),
#'                       v_ref = c(V1=30, V2=30, V3=0, V4=0),
#'                       lambda = c(0.5, 1, 1.5, 2), B = 10)
#' head(s$curve_simex)
#' }
#'
#' @export
simex_aft_spline <- function(
  data, x_var, sigma_w_sq,
  covariates = character(0), v_ref = NULL,
  df = 4, knots = NULL,
  outcome_var = "T_obs", status_var = "delta",
  lambda = c(0.5, 1, 1.5, 2), B = 50,
  dist = "lognormal", x_grid = NULL,
  verbose = FALSE, maxiter = 200
) {
  if (is.null(x_grid)) {
    x_grid <- seq(min(data[[x_var]]), max(data[[x_var]]), length.out = 100)
  }
  control_b <- survreg.control(maxiter = maxiter)
  n_fail <- 0L

  curve_at_lambda <- function(lam) {
    boots <- replicate(B, {
      data_b <- data
      data_b[[x_var]] <- data_b[[x_var]] +
        sqrt(lam * sigma_w_sq) * rnorm(nrow(data_b))
      converged <- TRUE
      fit_b <- withCallingHandlers(
        tryCatch(
          fit_aft_spline(data_b, x_var,
                         covariates = covariates,
                         df = df, knots = knots,
                         outcome_var = outcome_var, status_var = status_var,
                         dist = dist, control = control_b),
          error = function(e) NULL
        ),
        warning = function(w) {
          if (grepl("did not converge", w$message)) converged <<- FALSE
          invokeRestart("muffleWarning")
        }
      )
      if (is.null(fit_b) || !converged) {
        n_fail <<- n_fail + 1L
        return(rep(NA_real_, length(x_grid)))
      }
      lp <- predict_curve(fit_b, x_grid, v_ref, x_var)
      lp - lp[1]
    })
    rowMeans(boots, na.rm = TRUE)
  }

  fit_naive <- fit_aft_spline(data, x_var,
                              covariates = covariates,
                              df = df, knots = knots,
                              outcome_var = outcome_var, status_var = status_var,
                              dist = dist)
  lp_naive <- predict_curve(fit_naive, x_grid, v_ref, x_var)
  curve_naive <- lp_naive - lp_naive[1]

  pb <- if (verbose) txtProgressBar(min = 0, max = length(lambda), style = 3) else NULL
  curves_lambda <- sapply(seq_along(lambda), function(i) {
    y <- curve_at_lambda(lambda[i])
    if (!is.null(pb)) setTxtProgressBar(pb, i)
    y
  })
  if (!is.null(pb)) { close(pb); message("") }
  if (verbose && n_fail > 0L) {
    message(sprintf(
      "Dropped %d of %d perturbed fits (%.1f%%) for non-convergence.",
      n_fail, length(lambda) * B, 100 * n_fail / (length(lambda) * B)
    ))
  }

  lambda_full <- c(0, lambda)
  curves_full <- cbind(curve_naive, curves_lambda)

  curve_simex <- apply(curves_full, 1, function(yv) {
    df_extr <- data.frame(lam = lambda_full, y = yv)
    fit_extr <- lm(y ~ lam + I(lam^2), data = df_extr)
    as.numeric(predict(fit_extr, newdata = data.frame(lam = -1)))
  })

  list(x_grid = x_grid, curve_simex = curve_simex,
       curves_lambda = curves_full, lambda = lambda_full)
}
