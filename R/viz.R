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
#' @param ci_color Optional single colour for the confidence ribbon. If
#'   `NULL` (the default), uses `colors["SIMEX"]` when present (because the
#'   ribbon is typically a SIMEX bootstrap interval) and otherwise falls
#'   back to the first entry of `colors`.
#' @param time_ratio Logical; if `TRUE`, exponentiate the curves, truth, and
#'   CI ribbon so the y-axis reads as a time ratio relative to the reference
#'   anchor, and draw a horizontal reference line at `y = 1`. Default `FALSE`.
#' @param in_support Optional logical vector aligned to `x_grid` (e.g. from a
#'   [simex_aft_spline()] fit with `support_probs` set). Grid points where it
#'   is `FALSE` are drawn dashed with a lighter ribbon and a caption note,
#'   marking spline extrapolation beyond the exposure support.
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
                        colors = method_colors(),
                        ci_color = NULL,
                        time_ratio = FALSE,
                        in_support = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_curves(). ",
         "Install it via install.packages('ggplot2').")
  }

  # Drop entirely-NA series (a failed SIMEX/bootstrap curve) with a clear
  # message, rather than letting ggplot draw an empty panel and emit the
  # cryptic "no non-missing arguments to max" warning at render time.
  na_fit <- vapply(fits, function(f) all(!is.finite(f)), logical(1))
  if (any(na_fit)) {
    warning(sprintf(
      "Curve(s) %s are entirely NA and were dropped from the plot; the fit or bootstrap likely failed (check R_effective and any failure warnings).",
      paste(names(fits)[na_fit], collapse = ", ")), call. = FALSE)
    fits <- fits[!na_fit]
  }
  if (length(fits) == 0L) {
    stop("Nothing to plot: every curve is entirely NA. Inspect the SIMEX or ",
         "bootstrap output (e.g. 'R_effective' and the failure warnings) ",
         "before plotting.", call. = FALSE)
  }
  if (!is.null(ci) &&
      (all(!is.finite(ci$lower)) || all(!is.finite(ci$upper)))) {
    warning("Confidence interval is entirely NA; omitting the ribbon.",
            call. = FALSE)
    ci <- NULL
  }

  if (isTRUE(time_ratio)) {
    fits <- lapply(fits, exp)
    if (!is.null(truth)) truth <- exp(truth)
    if (!is.null(ci)) { ci$lower <- exp(ci$lower); ci$upper <- exp(ci$upper) }
    if (identical(y_lab, "Centred linear predictor (log time ratio)")) {
      y_lab <- "Time ratio (vs reference)"
    }
  }

  # Segments outside the exposure support are drawn dashed. `bridged` extends
  # the out-of-support mask by one neighbour on each side so the dashed and
  # solid portions share a vertex and the line has no visible gap.
  has_support <- !is.null(in_support) && any(!in_support)
  bridged <- NULL
  if (has_support) {
    n <- length(in_support)
    out <- !in_support
    bridged <- out
    bridged[-n] <- bridged[-n] | out[-1]
    bridged[-1] <- bridged[-1] | out[-n]
  }

  df_long <- do.call(rbind, lapply(names(fits), function(nm) {
    d <- data.frame(x = x_grid, y = fits[[nm]], method = nm)
    if (has_support) { d$in_support <- in_support; d$bridged <- bridged }
    d
  }))

  p <- ggplot2::ggplot()

  if (!is.null(ci)) {
    if (is.null(ci_color)) {
      ci_color <- if ("SIMEX" %in% names(colors)) unname(colors["SIMEX"]) else unname(colors[1])
    }
    df_ci <- data.frame(x = x_grid, lower = ci$lower, upper = ci$upper)
    if (has_support) {
      df_ci$in_support <- in_support
      df_ci$bridged    <- bridged
      p <- p + ggplot2::geom_ribbon(
        data = df_ci[df_ci$in_support, , drop = FALSE],
        ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
        fill = ci_color, alpha = 0.22)
      p <- p + ggplot2::geom_ribbon(
        data = df_ci[df_ci$bridged, , drop = FALSE],
        ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
        fill = ci_color, alpha = 0.08)
    } else {
      p <- p + ggplot2::geom_ribbon(
        data = df_ci,
        ggplot2::aes(x = .data$x, ymin = .data$lower, ymax = .data$upper),
        fill = ci_color, alpha = 0.22)
    }
  }

  if (!is.null(truth)) {
    df_truth <- data.frame(x = x_grid, y = truth)
    p <- p + ggplot2::geom_line(
      data = df_truth,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "black", linewidth = 1.0, linetype = "dashed"
    )
  }

  if (has_support) {
    p <- p +
      ggplot2::geom_line(
        data = df_long[df_long$in_support, , drop = FALSE],
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$method),
        linewidth = 0.9) +
      ggplot2::geom_line(
        data = df_long[df_long$bridged, , drop = FALSE],
        ggplot2::aes(x = .data$x, y = .data$y, color = .data$method,
                     group = .data$method),
        linewidth = 0.9, linetype = "dashed")
  } else {
    p <- p + ggplot2::geom_line(
      data = df_long,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$method),
      linewidth = 0.9)
  }

  if (isTRUE(time_ratio)) {
    p <- p + ggplot2::geom_hline(yintercept = 1, color = "grey50",
                                 linewidth = 0.4)
  }

  p <- p +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(x = x_lab, y = y_lab, title = title, color = "Method") +
    ggplot2::theme_bw(base_size = 12)

  if (has_support) {
    p <- p + ggplot2::labs(
      caption = "Dashed segments fall outside the exposure support (spline extrapolation, not fit).")
  }
  p
}
