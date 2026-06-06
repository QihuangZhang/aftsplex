test_that("plot_curves renders time-ratio and out-of-support variants", {
  skip_if_not_installed("ggplot2")
  x <- seq(6, 14, length.out = 40)
  truth <- f_true(x) - f_true(x[1])
  fits <- list(SIMEX = truth, Naive = 0.8 * truth)
  ci <- list(lower = truth - 0.1, upper = truth + 0.1)

  expect_s3_class(plot_curves(x, truth, fits, ci), "ggplot")
  expect_s3_class(plot_curves(x, truth, fits, ci, time_ratio = TRUE), "ggplot")

  in_support <- x >= 7 & x <= 13
  p <- plot_curves(x, truth, fits, ci, in_support = in_support)
  expect_s3_class(p, "ggplot")
  # out-of-support segments add extra layers (split solid/dashed + 2 ribbons)
  expect_gt(length(p$layers), length(plot_curves(x, truth, fits, ci)$layers))
})
