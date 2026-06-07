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
  if (length(sigma_u_sq) < 1L || !all(is.finite(sigma_u_sq)) ||
        any(sigma_u_sq <= 0)) {
    stop("'sigma_u_sq' must be a vector of finite positive variances.",
         call. = FALSE)
  }
  if (!is.finite(rho) || rho < -1 || rho > 1) {
    stop("'rho' must be a correlation in [-1, 1].", call. = FALSE)
  }
  sds <- sqrt(sigma_u_sq)
  R <- matrix(rho, length(sds), length(sds))
  diag(R) <- 1
  if (any(eigen(R, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("The implied correlation matrix is not positive definite for rho = ",
         rho, " with ", length(sds), " surrogates.", call. = FALSE)
  }
  outer(sds, sds) * R
}

#' Build a measurement-error covariance targeting a given reliability
#'
#' Returns a `Sigma_u` (via [build_sigma_u()]) whose surrogate-measurement
#' error variances are scaled so the **mean per-surrogate reliability** equals
#' `reliability`. Surrogate reliability is the intraclass-style ratio
#' `var_x / (var_x + var_U_j)`; the relative surrogate qualities in `ratios`
#' are preserved (so the surrogates stay unequal for the per-surrogate ladder),
#' and a single multiplicative scale is solved by [stats::uniroot()] to hit the
#' target average. Use it to sweep measurement quality in a simulation, e.g.
#' `reliability` in `c(0.5, 0.65, 0.8)`.
#'
#' @param reliability Target mean per-surrogate reliability in `(0, 1)`.
#' @param var_x Variance of the latent exposure `X` (e.g. `x_sd^2` passed to
#'   [generate_aft_data()]).
#' @param ratios Relative marginal error variances across the surrogates;
#'   their scale is solved for, only their ratios matter. Default
#'   `c(2.0, 2.5, 3.0)` matches [build_sigma_u()].
#' @param rho Common pairwise error correlation, passed to [build_sigma_u()].
#'
#' @return A `J x J` covariance matrix suitable for the `Sigma_u` argument of
#'   [generate_aft_data()].
#'
#' @examples
#' S <- sigma_u_for_reliability(0.65, var_x = 75^2)
#' # mean per-surrogate reliability is ~0.65:
#' mean(75^2 / (75^2 + diag(S)))
#'
#' @export
sigma_u_for_reliability <- function(reliability, var_x,
                                    ratios = c(2.0, 2.5, 3.0), rho = 0.4) {
  if (reliability <= 0 || reliability >= 1) {
    stop("'reliability' must be in (0, 1).", call. = FALSE)
  }
  if (length(var_x) != 1L || !is.finite(var_x) || var_x <= 0) {
    stop("'var_x' must be a single finite positive exposure variance.",
         call. = FALSE)
  }
  if (length(ratios) < 1L || !all(is.finite(ratios)) || any(ratios <= 0)) {
    stop("'ratios' must be a vector of finite positive relative variances.",
         call. = FALSE)
  }
  mean_rel <- function(scale) mean(var_x / (var_x + scale * ratios))
  # mean_rel is 1 at scale -> 0 and 0 at scale -> Inf, so a root exists.
  upper <- 1
  while (mean_rel(upper) > reliability) upper <- upper * 2
  scale <- stats::uniroot(function(s) mean_rel(s) - reliability,
                          interval = c(0, upper))$root
  build_sigma_u(scale * ratios, rho = rho)
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
#'   marginal variances `(2.0, 2.5, 3.0)` and `rho = 0.4`). Use
#'   [sigma_u_for_reliability()] to target a given reliability.
#' @param x_mean,x_sd Mean and standard deviation of the latent exposure
#'   `X_true` (shared by the main and validation samples). Defaults `10` and
#'   `sqrt(5)` reproduce the original arbitrary-unit scale; set e.g.
#'   `x_mean = 300`, `x_sd = 75` for a daily-minutes scale.
#' @param f_args Named list of extra arguments forwarded to [f_true()] for the
#'   dose-response (e.g. `list(x_min = 150, alpha = 0.01)` on a minutes scale).
#'   Default `list()` uses the [f_true()] defaults.
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
  x_mean     = 10,
  x_sd       = sqrt(5),
  f_args     = list(),
  mu         = 4,
  sigma_T    = 0.8,
  gamma      = c(0.10, -0.05, 0.30, -0.20),
  cens_rate  = 0.0005,
  seed       = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  if (length(x_sd) != 1L || !is.finite(x_sd) || x_sd < 0) {
    stop("'x_sd' must be a single finite non-negative standard deviation.",
         call. = FALSE)
  }
  J <- nrow(Sigma_u)
  f_eval <- function(x) do.call(f_true, c(list(x), f_args))

  X_true <- rnorm(n, mean = x_mean, sd = x_sd)
  V <- cbind(
    V1 = rnorm(n, 30, 5),
    V2 = rnorm(n, 30, 5),
    V3 = rbinom(n, 1, 0.7),
    V4 = rbinom(n, 1, 0.8)
  )

  log_T <- mu + f_eval(X_true) + V %*% gamma + sigma_T * rnorm(n)
  T_event <- exp(log_T)
  C_time  <- rexp(n, rate = cens_rate)
  T_obs   <- pmin(T_event, C_time)
  delta   <- as.integer(T_event <= C_time)

  U <- MASS::mvrnorm(n, mu = rep(0, J), Sigma = Sigma_u)
  W <- X_true + U
  colnames(W) <- paste0("W", seq_len(J))

  X_val <- rnorm(n_val, mean = x_mean, sd = x_sd)
  U_val <- MASS::mvrnorm(n_val, mu = rep(0, J), Sigma = Sigma_u)
  W_val <- X_val + U_val
  colnames(W_val) <- paste0("W", seq_len(J))

  list(
    survival = data.frame(T_obs = T_obs, delta = delta, X_true = X_true, W, V),
    validation = data.frame(X_true = X_val, W_val),
    truth = list(
      mu = mu, sigma_T = sigma_T, gamma = gamma,
      Sigma_u = Sigma_u, x_mean = x_mean, x_sd = x_sd,
      f_true = f_true, f_args = f_args
    )
  )
}
