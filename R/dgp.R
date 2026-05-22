#' True dose-response used by the simulation DGP
#'
#' Saturating exponential `beta * (1 - exp(-alpha * pmax(x - x_min, 0)))`.
#' This is the function the simulation tries to recover; it serves as the
#' ground truth in [generate_aft_data()] and is exported so users can
#' overlay it on plots.
#'
#' @param x Numeric vector of exposure values.
#' @param beta Saturation amplitude (asymptote as `x` grows large).
#' @param alpha Curvature parameter (larger = steeper initial rise).
#' @param x_min Threshold below which the response is flat at zero.
#'
#' @return Numeric vector the same length as `x`.
#'
#' @examples
#' x <- seq(0, 20, length.out = 100)
#' plot(x, f_true(x), type = "l")
#'
#' @export
f_true <- function(x, beta = 0.6, alpha = 0.4, x_min = 6) {
  beta * (1 - exp(-alpha * pmax(x - x_min, 0)))
}

#' Build a measurement-error covariance from marginal variances and a correlation
#'
#' Assembles a `J x J` covariance with the supplied marginal variances on the
#' diagonal and a common pairwise correlation `rho` off the diagonal. Used to
#' configure the multivariate surrogate noise in [generate_aft_data()].
#'
#' @param sigma_u_sq Numeric vector of marginal variances (length `J`).
#' @param rho Common pairwise correlation in `[-1/(J-1), 1]`.
#'
#' @return A `J x J` numeric matrix.
#'
#' @examples
#' build_sigma_u(c(2.0, 2.5, 3.0), rho = 0.4)
#'
#' @export
build_sigma_u <- function(sigma_u_sq = c(2.0, 2.5, 3.0), rho = 0.4) {
  sds <- sqrt(sigma_u_sq)
  R <- matrix(rho, length(sds), length(sds))
  diag(R) <- 1
  outer(sds, sds) * R
}

#' Generate one Monte Carlo replicate of the AFT-spline-SIMEX scenario
#'
#' Produces a main-study survival sample and an external validation sample
#' from the same data-generating process described in the manuscript
#' appendix. The latent log-exposure is normal, the dose-response is the
#' saturating curve [f_true()], the survival times follow a log-normal AFT
#' model with four confounders (V1-V4), and the surrogate measurements are
#' multivariate-normal around the latent truth with covariance `Sigma_u`.
#'
#' The returned `survival` data frame is in the shape expected by
#' [fit_aft_spline()] and [two_stage_bootstrap()] with their defaults
#' (`outcome_var = "T_obs"`, `status_var = "delta"`,
#' `covariates = c("V1", "V2", "V3", "V4")`,
#' `surrogate_cols = NULL` -> auto-detect via [extract_surrogate_cols()]).
#'
#' @param n Main-study sample size.
#' @param n_val Validation sample size.
#' @param Sigma_u `J x J` measurement-error covariance (default 3x3 with
#'   marginal variances `(2.0, 2.5, 3.0)` and `rho = 0.4`).
#' @param mu AFT intercept.
#' @param sigma_T AFT log-scale standard deviation.
#' @param gamma Numeric vector of length 4: confounder coefficients in the
#'   linear predictor (`V1, V2, V3, V4`).
#' @param cens_rate Rate of the exponential censoring distribution. Smaller
#'   values give heavier event experience. The default `0.0005` yields
#'   ~76% events under the other defaults.
#' @param seed Optional integer seed.
#'
#' @return A list with elements `survival` (main-study data frame),
#'   `validation` (validation data frame with paired truth + surrogates),
#'   and `truth` (a list of the parameter values used).
#'
#' @examples
#' sim <- generate_aft_data(n = 200, n_val = 100, seed = 1)
#' names(sim)
#' head(sim$survival)
#'
#' @export
generate_aft_data <- function(
  n          = 2000,
  n_val      = 500,
  Sigma_u    = build_sigma_u(),
  mu         = 4,
  sigma_T    = 0.8,
  gamma      = c(0.10, -0.05, 0.30, -0.20),
  cens_rate  = 0.0005,
  seed       = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  J <- nrow(Sigma_u)

  X_true <- rnorm(n, mean = 10, sd = sqrt(5))
  V <- cbind(
    V1 = rnorm(n, 30, 5),
    V2 = rnorm(n, 30, 5),
    V3 = rbinom(n, 1, 0.7),
    V4 = rbinom(n, 1, 0.8)
  )

  log_T <- mu + f_true(X_true) + V %*% gamma + sigma_T * rnorm(n)
  T_event <- exp(log_T)
  C_time  <- rexp(n, rate = cens_rate)
  T_obs   <- pmin(T_event, C_time)
  delta   <- as.integer(T_event <= C_time)

  U <- MASS::mvrnorm(n, mu = rep(0, J), Sigma = Sigma_u)
  W <- X_true + U
  colnames(W) <- paste0("W", seq_len(J))

  X_val <- rnorm(n_val, mean = 10, sd = sqrt(5))
  U_val <- MASS::mvrnorm(n_val, mu = rep(0, J), Sigma = Sigma_u)
  W_val <- X_val + U_val
  colnames(W_val) <- paste0("W", seq_len(J))

  list(
    survival = data.frame(T_obs = T_obs, delta = delta, X_true = X_true, W, V),
    validation = data.frame(X_true = X_val, W_val),
    truth = list(
      mu = mu, sigma_T = sigma_T, gamma = gamma,
      Sigma_u = Sigma_u, f_true = f_true
    )
  )
}
