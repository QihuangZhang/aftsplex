#' Fit the Phase-1 multivariate ME calibration model
#'
#' Fits per-equation OLS of each surrogate column on the validation truth:
#' `W_j = alpha_0j + alpha_1j * X + e_j`, with `e ~ N_J(0, Sigma_e)`.
#' Under shared regressors this is equivalent to seemingly-unrelated
#' regression (SUR) and admits closed-form ML.
#'
#' @param validation Validation data frame containing the truth column named
#'   in `x_var` and the surrogate columns matched by `surrogate_pattern`.
#' @param x_var Name of the truth column. Default `"X_true"`.
#' @param surrogate_pattern Regex used to identify surrogate columns; see
#'   [extract_surrogate_cols()]. Default `"^W[0-9]+$"`.
#'
#' @return A list with elements
#'   * `alpha0` numeric vector of intercepts (length `J`),
#'   * `alpha1` numeric vector of slopes (length `J`),
#'   * `Sigma_e` `J x J` residual covariance,
#'   * `W_cols` character vector of surrogate column names used,
#'   * `n_val` number of validation rows.
#'
#' @examples
#' sim <- generate_aft_data(n = 100, n_val = 200, seed = 1)
#' fit_me_calibration(sim$validation)
#'
#' @export
fit_me_calibration <- function(validation, x_var = "X_true",
                               surrogate_pattern = "^W[0-9]+$") {
  W_cols <- extract_surrogate_cols(validation, surrogate_pattern)
  J <- length(W_cols)
  X <- validation[[x_var]]
  if (is.null(X)) {
    stop("Truth column '", x_var, "' not found in 'validation'.", call. = FALSE)
  }
  if (nrow(validation) < 3L) {
    stop("'validation' has ", nrow(validation), " row(s); at least 3 are needed ",
         "to identify the calibration slopes and a residual variance.",
         call. = FALSE)
  }
  if (!all(is.finite(X)) || stats::var(X) < .Machine$double.eps) {
    stop("The validation truth column '", x_var, "' is constant or non-finite; ",
         "the calibration slopes are unidentified.", call. = FALSE)
  }
  alpha0 <- numeric(J)
  alpha1 <- numeric(J)
  resid_mat <- matrix(NA_real_, nrow(validation), J)
  for (j in seq_len(J)) {
    fit_j <- lm(validation[[W_cols[j]]] ~ X)
    alpha0[j] <- coef(fit_j)[1]
    alpha1[j] <- coef(fit_j)[2]
    resid_mat[, j] <- residuals(fit_j)
  }
  Sigma_e <- crossprod(resid_mat) / (nrow(validation) - 2)
  list(alpha0 = alpha0, alpha1 = alpha1, Sigma_e = Sigma_e,
       W_cols = W_cols, n_val = nrow(validation))
}

#' GLS-combine surrogate columns into a single calibrated exposure
#'
#' Applies the linear back-transform `W_tilde = (W - alpha_0) / alpha_1` and
#' the GLS combiner `omega = Sigma_tilde^-1 1 / (1' Sigma_tilde^-1 1)` to
#' produce a single calibrated exposure with conditional variance
#' `sigma_w_sq = (1' Sigma_tilde^-1 1)^-1`.
#'
#' @param W Numeric matrix of surrogate columns, columns ordered as in
#'   `calibration$W_cols`.
#' @param calibration Output of [fit_me_calibration()].
#'
#' @return A list with elements
#'   * `W_bar` numeric vector of length `nrow(W)`,
#'   * `sigma_w_sq` scalar conditional variance,
#'   * `omega` numeric vector of GLS weights (length `J`),
#'   * `Sigma_tilde` `J x J` back-transformed error covariance.
#'
#' @examples
#' sim <- generate_aft_data(n = 100, n_val = 200, seed = 1)
#' cal <- fit_me_calibration(sim$validation)
#' W   <- as.matrix(sim$survival[, cal$W_cols])
#' out <- gls_combine(W, cal)
#' head(out$W_bar)
#'
#' @export
gls_combine <- function(W, calibration) {
  alpha0 <- calibration$alpha0
  alpha1 <- calibration$alpha1
  Sigma_e <- calibration$Sigma_e
  J <- length(alpha1)

  if (!all(is.finite(W))) {
    stop("'W' contains non-finite values (NA/NaN/Inf); cannot GLS-combine the ",
         "surrogates.", call. = FALSE)
  }
  if (any(abs(alpha1) < 1e-8) || !all(is.finite(alpha1))) {
    stop("A calibration slope ('alpha1') is ~0 or non-finite, so the surrogate ",
         "carries no information about the exposure and cannot be ",
         "back-transformed. Drop the uninformative surrogate.", call. = FALSE)
  }

  W_tilde <- sweep(sweep(W, 2, alpha0, "-"), 2, alpha1, "/")
  D_inv <- diag(1 / alpha1, nrow = J)
  Sigma_tilde <- D_inv %*% Sigma_e %*% D_inv
  Sigma_tilde_inv <- tryCatch(
    solve(Sigma_tilde),
    error = function(e) {
      stop("The back-transformed surrogate covariance is singular; this is ",
           "usually collinear or constant surrogate columns. Drop a redundant ",
           "surrogate.", call. = FALSE)
    }
  )
  ones <- rep(1, J)
  denom <- as.numeric(t(ones) %*% Sigma_tilde_inv %*% ones)
  omega <- as.vector(Sigma_tilde_inv %*% ones) / denom
  W_bar <- as.vector(W_tilde %*% omega)
  sigma_w_sq <- 1 / denom

  if (!is.finite(sigma_w_sq) || sigma_w_sq < 0) {
    stop("The combined conditional variance 'sigma_w_sq' is non-finite or ",
         "negative; the calibration is degenerate.", call. = FALSE)
  }

  list(W_bar = W_bar, sigma_w_sq = sigma_w_sq, omega = omega,
       Sigma_tilde = Sigma_tilde)
}
