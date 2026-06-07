#' Integrated squared error of a fitted curve against the truth
#'
#' Trapezoidal-rule approximation of `integral of (f_hat - f_true)^2 dx`
#' over the exposure grid. The truth function is centred at the lower grid
#' boundary to match the convention used elsewhere in the package
#' (`truth(x_grid) - truth(x_grid[1])`).
#'
#' @param x_grid Numeric vector of exposure values.
#' @param truth_fn A function such as [f_true()] returning the (uncentred)
#'   true response at `x_grid`.
#' @param fits A named list of numeric vectors, each of length
#'   `length(x_grid)`, holding fitted centred curves.
#'
#' @return A list with elements `truth` (the centred truth at `x_grid`) and
#'   `ise` (named numeric vector of ISEs, one entry per element of `fits`).
#'
#' @examples
#' x_grid <- seq(6, 14, length.out = 50)
#' fits <- list(naive = numeric(length(x_grid)),
#'              oracle = f_true(x_grid) - f_true(x_grid[1]))
#' compute_ise(x_grid, f_true, fits)
#'
#' @export
compute_ise <- function(x_grid, truth_fn, fits) {
  if (length(x_grid) < 2L) {
    stop("'x_grid' must have at least 2 points to integrate.", call. = FALSE)
  }
  if (is.unsorted(x_grid)) {
    stop("'x_grid' must be sorted in increasing order; otherwise the ",
         "trapezoidal integral is undefined (and can come out negative).",
         call. = FALSE)
  }
  bad <- vapply(fits, function(f) length(f) != length(x_grid), logical(1))
  if (any(bad)) {
    stop("Each element of 'fits' must have length(x_grid) = ", length(x_grid),
         "; offending element(s): ",
         paste(names(fits)[bad], collapse = ", "), ".", call. = FALSE)
  }
  truth <- truth_fn(x_grid) - truth_fn(x_grid[1])
  ise <- sapply(fits, function(f) {
    delta_x <- diff(x_grid)
    sq_err <- ((f - truth)[-1] + (f - truth)[-length(f)])^2 / 4
    sum(sq_err * delta_x)
  })
  list(truth = truth, ise = ise)
}
