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
#' @param x_ref Optional scalar reference exposure at which to anchor the
#'   centred curves (so the dose-response is read as relative to `x_ref`,
#'   e.g. a clinically meaningful value). `NULL` (default) anchors at the
#'   lower grid boundary `x_grid[1]`, preserving the original behaviour.
#' @param support_probs Optional length-2 numeric of probabilities defining
#'   the trustworthy exposure support as quantiles of `data[[x_var]]` (e.g.
#'   `c(0.05, 0.95)`). When set, grid points outside the support are flagged
#'   in the returned `in_support` vector and a one-shot `warning()` notes
#'   that they are spline extrapolations. `NULL` (default) disables the check.
#' @param verbose Print progress bar and non-convergence count?
#' @param maxiter `survreg` iteration cap (passed via
#'   [survival::survreg.control()]).
#'
#' @return An object of class `"simex_aft_spline"`: a list with elements
#'   * `x_grid` the exposure grid,
#'   * `curve_simex` the centred SIMEX-corrected linear predictor,
#'   * `curve_naive` the centred naive linear predictor (lambda = 0),
#'   * `curves_lambda` matrix of centred lambda-fits (columns = lambda, including 0),
#'   * `lambda` the lambda values used (with 0 in column 1),
#'   * `in_support` logical vector aligned to `x_grid` (or `NULL` if
#'     `support_probs` was not supplied),
#'   * `x_ref` the reference anchor used (or `NULL`),
#'   * `call`, `n`, `df`, `dist`, `covariates`, `B`, `sigma_w_sq`, `n_fail`
#'     configuration and diagnostics used by [summary.simex_aft_spline()]
#'     and [plot.simex_aft_spline()].
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
#' summary(s)
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
  x_ref = NULL, support_probs = NULL,
  verbose = FALSE, maxiter = 200
) {
  call <- match.call()
  if (is.null(data[[x_var]])) {
    stop("Exposure column '", x_var, "' not found in 'data'.", call. = FALSE)
  }
  if (!is.null(x_grid) && length(x_grid) < 2L) {
    stop("'x_grid' must have at least 2 points.", call. = FALSE)
  }
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

  fit_naive <- tryCatch(
    fit_aft_spline(data, x_var,
                   covariates = covariates,
                   df = df, knots = knots,
                   outcome_var = outcome_var, status_var = status_var,
                   dist = dist),
    error = function(e) {
      stop("SIMEX naive AFT fit failed: ", conditionMessage(e),
           ". This usually means too few events or a rank-deficient spline ",
           "basis; reduce 'df' or check the data.", call. = FALSE)
    }
  )
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
  colnames(curves_full) <- paste0("lam=", format(lambda_full, digits = 2))

  # Pointwise quadratic extrapolation to lambda = -1. Grid points whose
  # lambda-curve has fewer than three finite values (e.g. perturbed fits that
  # failed to converge) cannot support the quadratic and are returned as NA
  # rather than triggering a rank-deficient fit or an "0 (non-NA) cases"
  # error; the count is reported via a single warning.
  n_unextrap <- 0L
  curve_simex <- apply(curves_full, 1, function(yv) {
    ok <- is.finite(yv)
    if (sum(ok) < 3L) {
      n_unextrap <<- n_unextrap + 1L
      return(NA_real_)
    }
    df_extr <- data.frame(lam = lambda_full[ok], y = yv[ok])
    fit_extr <- lm(y ~ lam + I(lam^2), data = df_extr)
    as.numeric(predict(fit_extr, newdata = data.frame(lam = -1)))
  })
  if (n_unextrap > 0L) {
    warning(sprintf(
      paste0("SIMEX extrapolation undefined at %d of %d grid point(s) ",
             "(too few converged perturbed fits); returned NA there."),
      n_unextrap, length(x_grid)), call. = FALSE)
  }

  if (!is.null(x_ref)) {
    curve_simex <- recenter_curve(curve_simex, x_grid, x_ref)
    curve_naive <- recenter_curve(curve_naive, x_grid, x_ref)
    cn <- colnames(curves_full)
    curves_full <- apply(curves_full, 2, recenter_curve,
                         x_grid = x_grid, x_ref = x_ref)
    colnames(curves_full) <- cn
  }

  in_support <- support_flag(x_grid, data[[x_var]], support_probs)

  out <- list(
    x_grid       = x_grid,
    curve_simex  = curve_simex,
    curve_naive  = curve_naive,
    curves_lambda = curves_full,
    lambda       = lambda_full,
    in_support   = in_support,
    x_ref        = x_ref,
    call         = call,
    n            = nrow(data),
    df           = df,
    dist         = dist,
    covariates   = covariates,
    B            = B,
    sigma_w_sq   = sigma_w_sq,
    n_fail       = n_fail
  )
  class(out) <- "simex_aft_spline"
  out
}

#' @export
print.simex_aft_spline <- function(x, ...) {
  cat("SIMEX-corrected AFT-spline dose-response\n\n")
  cat("Call:\n  "); print(x$call); cat("\n")
  cat(sprintf("  x_grid: [%.3f, %.3f] (%d points)   n_obs: %d   df: %d\n",
              min(x$x_grid), max(x$x_grid), length(x$x_grid), x$n, x$df))
  cat(sprintf("  sigma_w_sq: %.4f   lambda: %s   B: %d\n",
              x$sigma_w_sq,
              paste(format(x$lambda[x$lambda > 0], digits = 2), collapse = ","),
              x$B))
  cat("\nUse summary() for a quantile table, plot() to draw the curve.\n")
  invisible(x)
}

#' Summarise a SIMEX-corrected AFT-spline fit
#'
#' Returns a tabular summary of the SIMEX point estimate together with the
#' configuration that produced it. The headline table reports the centred
#' Naive and SIMEX dose-response (and their difference, the SIMEX
#' correction) at quantiles of the exposure grid.
#'
#' @param object A [simex_aft_spline()] fit.
#' @param probs Quantiles of `x_grid` at which to report the centred curves.
#'   Default `c(0.05, 0.25, 0.50, 0.75, 0.95)`.
#' @param ... Unused.
#'
#' @return An object of class `"summary.simex_aft_spline"` with a `table`
#'   element and the configuration fields. Has a `print` method.
#'
#' @export
summary.simex_aft_spline <- function(object,
                                     probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
                                     ...) {
  qs <- as.numeric(quantile(object$x_grid, probs))
  approx_curve <- function(yv) stats::approx(object$x_grid, yv, xout = qs)$y
  naive <- object$curve_naive
  simex <- object$curve_simex
  tab <- rbind(
    x          = qs,
    Naive      = approx_curve(naive),
    SIMEX      = approx_curve(simex),
    Correction = approx_curve(simex - naive)
  )
  colnames(tab) <- sprintf("%g%%", 100 * probs)

  out <- list(
    call       = object$call,
    n          = object$n,
    df         = object$df,
    dist       = object$dist,
    covariates = object$covariates,
    lambda     = object$lambda,
    B          = object$B,
    sigma_w_sq = object$sigma_w_sq,
    n_fail     = object$n_fail,
    table      = tab
  )
  class(out) <- "summary.simex_aft_spline"
  out
}

#' @export
print.summary.simex_aft_spline <- function(x, digits = 3, ...) {
  cat("SIMEX-corrected AFT-spline dose-response\n\n")
  cat("Call:\n  "); print(x$call); cat("\n")
  cat("Configuration:\n")
  cov_str <- if (length(x$covariates) == 0L) "(none)" else paste(x$covariates, collapse = ", ")
  cat(sprintf("  n observations:  %-18d Spline df:    %d\n", x$n, x$df))
  cat(sprintf("  Covariates:      %-18s Distribution: %s\n", cov_str, x$dist))
  cat(sprintf("  Lambda grid:     %s\n",
              paste(format(x$lambda, digits = 2), collapse = ", ")))
  cat(sprintf("  Inner B:         %-18d sigma_w_sq:   %.4f\n",
              x$B, x$sigma_w_sq))
  cat("\nCentred dose-response (anchored at min(x_grid)):\n")
  print(round(x$table, digits))
  total <- x$B * sum(x$lambda > 0)
  if (!is.null(x$n_fail) && x$n_fail > 0L) {
    cat(sprintf("\nNotes: %d of %d perturbed fits dropped for non-convergence.\n",
                x$n_fail, total))
  }
  invisible(x)
}

#' Plot a SIMEX-corrected AFT-spline fit
#'
#' Convenience S3 method that draws the SIMEX-corrected curve via
#' [plot_curves()]. By default the naive curve (`lambda = 0`) is overlaid
#' as a smoothing-bias diagnostic. Pass `truth` to add the dashed truth
#' line, and `ci` to add a percentile ribbon (typically from
#' [two_stage_bootstrap()]).
#'
#' @param x A [simex_aft_spline()] fit.
#' @param truth Optional numeric vector (length `length(x$x_grid)`) of the
#'   centred truth curve to overlay as a dashed black line.
#' @param ci Optional list with elements `lower` and `upper`.
#' @param show_naive Logical; overlay the naive (`lambda = 0`) curve?
#'   Default `TRUE`.
#' @param time_ratio Logical; plot the exponentiated curve (time ratio
#'   relative to the reference anchor) instead of the centred linear
#'   predictor? Default `FALSE`. Passed to [plot_curves()].
#' @param title Plot title.
#' @param ... Passed through to [plot_curves()].
#'
#' @return A `ggplot` object.
#'
#' @export
plot.simex_aft_spline <- function(x, truth = NULL, ci = NULL,
                                  show_naive = TRUE, time_ratio = FALSE,
                                  title = "SIMEX-corrected AFT-spline dose-response",
                                  ...) {
  fits <- list(SIMEX = x$curve_simex)
  if (isTRUE(show_naive)) fits$Naive <- x$curve_naive
  plot_curves(x$x_grid, truth = truth, fits = fits, ci = ci,
              title = title, time_ratio = time_ratio,
              in_support = x$in_support, ...)
}
