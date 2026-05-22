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
  expect_named(out, c("x_grid", "curve_simex", "curves_lambda", "lambda"))
  expect_length(out$curve_simex, 25L)
  expect_equal(dim(out$curves_lambda), c(25L, 4L))   # 1 naive + 3 lambdas
  expect_equal(out$lambda[1], 0)
  expect_true(all(is.finite(out$curve_simex)))
})
