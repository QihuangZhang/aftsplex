test_that("simex_aft_spline returns expected list shape", {
  sim <- generate_aft_data(n = 300, n_val = 100, seed = 1)
  x_grid <- seq(7, 13, length.out = 25)
  set.seed(123)
  out <- simex_aft_spline(
    sim$survival, x_var = "X_true",
    sigma_w_sq = 1.0,
    covariates = c("V1","V2","V3","V4"),
    v_ref = c(V1 = 30, V2 = 30, V3 = 0, V4 = 0),
    lambda = c(0.5, 1, 2), B = 3,
    x_grid = x_grid
  )
  expect_s3_class(out, "simex_aft_spline")
  expect_true(all(c("x_grid", "curve_simex", "curve_naive", "curves_lambda",
                    "lambda", "call", "n", "df", "dist", "covariates",
                    "B", "sigma_w_sq", "n_fail") %in% names(out)))
  expect_length(out$curve_simex, 25L)
  expect_length(out$curve_naive, 25L)
  expect_equal(dim(out$curves_lambda), c(25L, 4L))   # 1 naive + 3 lambdas
  expect_equal(out$lambda[1], 0)
  expect_true(all(is.finite(out$curve_simex)))
})

test_that("simex_aft_spline S3 methods dispatch", {
  sim <- generate_aft_data(n = 300, n_val = 100, seed = 1)
  set.seed(123)
  out <- simex_aft_spline(
    sim$survival, x_var = "X_true",
    sigma_w_sq = 1.0,
    covariates = c("V1","V2","V3","V4"),
    v_ref = c(V1 = 30, V2 = 30, V3 = 0, V4 = 0),
    lambda = c(0.5, 1), B = 2,
    x_grid = seq(7, 13, length.out = 10)
  )
  expect_output(print(out), "SIMEX")
  s <- summary(out)
  expect_s3_class(s, "summary.simex_aft_spline")
  expect_true(all(c("x", "Naive", "SIMEX", "Correction") %in% rownames(s$table)))
  expect_output(print(s), "Centred dose-response")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(out), "ggplot")
  expect_s3_class(plot(out, time_ratio = TRUE), "ggplot")
})

test_that("x_ref and support_probs behave in simex_aft_spline", {
  sim <- generate_aft_data(n = 300, n_val = 100, seed = 1)
  x_grid <- seq(2, 18, length.out = 25)
  set.seed(123)
  expect_warning(
    out <- simex_aft_spline(
      sim$survival, x_var = "X_true", sigma_w_sq = 1.0,
      covariates = c("V1","V2","V3","V4"),
      v_ref = c(V1 = 30, V2 = 30, V3 = 0, V4 = 0),
      lambda = c(0.5, 1), B = 3, x_grid = x_grid,
      x_ref = x_grid[13], support_probs = c(0.05, 0.95)
    ),
    "extrapolations"
  )
  expect_equal(unname(out$curve_simex[13]), 0, tolerance = 1e-10)
  expect_length(out$in_support, length(x_grid))
  expect_false(out$in_support[1])
})
