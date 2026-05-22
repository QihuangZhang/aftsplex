#' aftsplex: AFT + spline + SIMEX for nonlinear dose-response with multivariate ME
#'
#' Two-phase estimator combining (i) a multivariate measurement-error
#' calibration on an external validation sample with a generalised-least-
#' squares (GLS) combiner, and (ii) a natural-cubic-spline accelerated-
#' failure-time model on the calibrated exposure, corrected for residual
#' attenuation via SIMEX. Inference uses a two-stage nonparametric bootstrap
#' that jointly resamples the validation and main-study samples.
#'
#' @keywords internal
#' @importFrom stats lm coef residuals predict rnorm rbinom rexp qnorm
#'   quantile sd median as.formula
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom survival survreg Surv survreg.control
#' @importFrom splines ns
#' @importFrom MASS mvrnorm
"_PACKAGE"

utils::globalVariables(c("x", "y", "method", "lower", "upper", ".data"))
