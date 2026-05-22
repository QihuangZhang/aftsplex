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
