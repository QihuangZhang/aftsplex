#' Fit an AFT model with a natural-cubic-spline exposure
#'
#' Thin wrapper around [survival::survreg()] that builds a formula of the
#' form `Surv(outcome_var, status_var) ~ ns(x_var, df = ...) + covariates`
#' from the column names supplied. Used three ways in the simulation: as the
#' oracle (on the true exposure), as the naive estimator (on the GLS-combined
#' surrogate without ME correction), and as the inner workhorse of
#' [simex_aft_spline()].
#'
#' @param data A data frame containing the outcome, status, exposure, and
#'   covariate columns named below.
#' @param x_var Name of the exposure column.
#' @param covariates Character vector of additional covariate column names.
#'   Pass `character(0)` (the default) for an exposure-only model.
#' @param df Spline degrees of freedom (ignored if `knots` is supplied).
#' @param knots Optional numeric vector of interior knot locations.
#' @param outcome_var Name of the time-to-event column. Default `"T_obs"`.
#' @param status_var Name of the event indicator column. Default `"delta"`.
#' @param dist `survreg` distribution. Default `"lognormal"` matches the DGP
#'   in [generate_aft_data()].
#' @param control Optional [survival::survreg.control()] object.
#'
#' @return A `survreg` fitted-model object.
#'
#' @examples
#' sim <- generate_aft_data(n = 200, n_val = 100, seed = 1)
#' fit <- fit_aft_spline(sim$survival, x_var = "X_true",
#'                       covariates = c("V1","V2","V3","V4"))
#' coef(fit)
#'
#' @export
fit_aft_spline <- function(data, x_var,
                           covariates = character(0),
                           df = 4, knots = NULL,
                           outcome_var = "T_obs",
                           status_var  = "delta",
                           dist = "lognormal", control = NULL) {
  if (is.null(knots)) {
    basis_call <- substitute(ns(z, df = D), list(z = as.name(x_var), D = df))
  } else {
    basis_call <- substitute(ns(z, knots = K), list(z = as.name(x_var), K = knots))
  }
  rhs <- paste(c(deparse(basis_call), covariates), collapse = " + ")
  fml <- as.formula(paste0(
    "Surv(", outcome_var, ", ", status_var, ") ~ ", rhs
  ))
  if (is.null(control)) {
    survreg(fml, data = data, dist = dist)
  } else {
    survreg(fml, data = data, dist = dist, control = control)
  }
}

#' Predict the centred dose-response curve at a reference covariate setting
#'
#' Returns the linear predictor of a fitted [fit_aft_spline()] model on an
#' exposure grid, holding the covariates fixed at `v_ref`. Subtracting the
#' first element produces the centred curve used throughout the manuscript.
#'
#' @param fit A `survreg` object returned by [fit_aft_spline()].
#' @param x_grid Numeric vector of exposure values to predict at.
#' @param v_ref Named numeric vector of reference covariate values, with
#'   names matching the `covariates` used when fitting. Pass `NULL` if the
#'   model has no covariates.
#' @param x_var Name of the exposure column (must match `fit_aft_spline`).
#'
#' @return A numeric vector of length `length(x_grid)` (the uncentred linear
#'   predictor; subtract `out[1]` to centre at the lower grid boundary).
#'
#' @examples
#' sim <- generate_aft_data(n = 300, n_val = 100, seed = 1)
#' fit <- fit_aft_spline(sim$survival, x_var = "X_true",
#'                       covariates = c("V1","V2","V3","V4"))
#' x_grid <- seq(7, 13, length.out = 50)
#' lp <- predict_curve(fit, x_grid,
#'                     v_ref = c(V1=30, V2=30, V3=0, V4=0),
#'                     x_var = "X_true")
#'
#' @export
predict_curve <- function(fit, x_grid, v_ref = NULL, x_var) {
  if (is.null(v_ref) || length(v_ref) == 0L) {
    newdata <- data.frame(x = x_grid)
    names(newdata) <- x_var
  } else {
    newdata <- data.frame(matrix(rep(v_ref, length(x_grid)),
                                 nrow = length(x_grid), byrow = TRUE))
    colnames(newdata) <- names(v_ref)
    newdata[[x_var]] <- x_grid
  }
  predict(fit, newdata = newdata, type = "lp")
}
