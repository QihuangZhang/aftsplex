#' Identify surrogate-measurement columns in a data frame
#'
#' Returns the names of columns matching a regex pattern. The default
#' `"^W[0-9]+$"` matches `W1`, `W2`, ... and corresponds to the convention
#' used by [generate_aft_data()] and the manuscript notation.
#'
#' @param data A data frame.
#' @param pattern A regular expression. Defaults to `"^W[0-9]+$"`.
#'
#' @return A character vector of column names, in the order they appear in
#'   `data`. Throws an error if no columns match.
#'
#' @examples
#' sim <- generate_aft_data(n = 50, n_val = 20, seed = 1)
#' extract_surrogate_cols(sim$survival)
#'
#' @export
extract_surrogate_cols <- function(data, pattern = "^W[0-9]+$") {
  cols <- grep(pattern, names(data), value = TRUE)
  if (length(cols) == 0L) {
    stop("No surrogate columns matched pattern '", pattern,
         "' in data with columns: ",
         paste(names(data), collapse = ", "))
  }
  cols
}

## Re-anchor a centred curve so that it equals zero at `x_ref`. Centring is an
## additive per-curve constant, so this is exact regardless of where the curve
## was originally centred. `NULL` x_ref is a no-op (keeps the grid-1 anchor).
recenter_curve <- function(curve, x_grid, x_ref) {
  if (is.null(x_ref)) return(curve)
  offset <- stats::approx(x_grid, curve, xout = x_ref, rule = 2)$y
  curve - offset
}

## Logical "is this grid point inside the exposure support?" plus a one-shot
## extrapolation warning. `support_probs` NULL disables the check entirely.
support_flag <- function(x_grid, x_obs, support_probs) {
  if (is.null(support_probs)) return(NULL)
  bounds <- stats::quantile(x_obs, probs = support_probs, names = FALSE,
                            na.rm = TRUE)
  in_support <- x_grid >= bounds[1] & x_grid <= bounds[2]
  n_out <- sum(!in_support)
  if (n_out > 0L) {
    warning(sprintf(
      paste0("%d of %d x_grid points fall outside the [%.3g, %.3g] support ",
             "(the %g--%g%% percentiles of the exposure); these are spline ",
             "extrapolations, not fits."),
      n_out, length(x_grid), bounds[1], bounds[2],
      100 * support_probs[1], 100 * support_probs[2]
    ), call. = FALSE)
  }
  in_support
}
