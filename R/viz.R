#' Overlay fitted dose-response curves with the truth and an optional CI
#'
#' Builds a faceted ggplot2 figure comparing a named list of fitted curves
#' against the truth, with an optional bootstrap confidence ribbon. Uses
#' the [method_colors()] palette by default.
#'
#' Requires the suggested package `ggplot2` to be installed.
#'
#' @param x_grid Numeric exposure grid.
#' @param truth Numeric vector of truth values (already centred at the
#'   first grid point), the same length as `x_grid`. Pass `NULL` to omit
#'   the truth line (e.g. on real data where the truth is unknown).
#' @param fits Named list of numeric vectors of length `length(x_grid)`.
#' @param ci Optional list with elements `lower` and `upper`, each a
#'   numeric vector of length `length(x_grid)`.
#' @param title Plot title.
#' @param x_lab,y_lab Axis labels.
#' @param colors Optional named character vector of colours; defaults to
#'   [method_colors()] (matching `Naive`, `Oracle`, `SIMEX`).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- seq(6, 14, length.out = 50)
#'   truth <- f_true(x) - f_true(x[1])
#'   plot_curves(x, truth, list(Oracle = truth))
#' }
#'
#' @export
plot_curves <- function(x_grid, truth, fits, ci = NULL,
                        title = "AFT + spline dose-response",
                        x_lab = "Exposure",
                        y_lab = "Centred linear predictor (log time ratio)",
                        colors = method_colors()) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_curves(). ",
         "Install it via install.packages('ggplot2').")
  }
  df_long <- do.call(rbind, lapply(names(fits), function(nm) {
    data.frame(x = x_grid, y = fits[[nm]], method = nm)
  }))

  p <- ggplot2::ggplot()
  if (!is.null(truth)) {
    df_truth <- data.frame(x = x_grid, y = truth)
    p <- p + ggplot2::geom_line(
      data = df_truth,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black", linewidth = 1.0, linetype = "dashed"
    )
  }
  p <- p +
    ggplot2::geom_line(data = df_long,
                       ggplot2::aes(x = .data$x, y = .data$y, color = .data$method),
                       linewidth = 0.9) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = x_lab, y = y_lab, title = title, color = "Method") +
    ggplot2::theme_bw(base_size = 12)

  if (!is.null(ci)) {
    df_ci <- data.frame(x = x_grid, lower = ci$lower, upper = ci$upper)
    p <- p + ggplot2::geom_ribbon(
      data = df_ci,
      ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
      fill = unname(colors[1]), alpha = 0.18
    )
  }
  p
}
