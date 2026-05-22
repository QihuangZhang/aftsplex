#' StaGILL lab palette for this project
#'
#' Named hex colours used by [plot_curves()] and the package vignettes.
#' Project-specific palette for the AFT-spline-SIMEX manuscript (see
#' `feedback-lab-palette` memo): every plot in the package uses an
#' explicit named palette via `scale_*_manual()` rather than the ggplot2
#' default.
#'
#' @return A named character vector of hex colours.
#' @examples
#' stagill_palette()
#' @export
stagill_palette <- function() {
  c(
    Eggshell       = "#f4f1de",
    BurntPeach     = "#e07a5f",
    TwilightIndigo = "#3d405b",
    MutedTeal      = "#81b29a",
    ApricotCream   = "#f2cc8f"
  )
}

#' Semantic method-to-colour mapping for Naive / Oracle / SIMEX comparisons
#'
#' Maps the three estimator names used throughout the package to colours
#' drawn from [stagill_palette()]. Pass to `ggplot2::scale_color_manual()`.
#'
#' @return A named character vector with entries `Naive`, `Oracle`, `SIMEX`.
#' @examples
#' method_colors()
#' @export
method_colors <- function() {
  p <- stagill_palette()
  c(Naive  = unname(p["BurntPeach"]),
    Oracle = unname(p["TwilightIndigo"]),
    SIMEX  = unname(p["MutedTeal"]))
}
