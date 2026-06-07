test_that("f_true is zero below threshold and positive above", {
  expect_equal(f_true(0), 0)
  expect_equal(f_true(6), 0)
  expect_gt(f_true(10), 0)
  expect_lt(f_true(100), 0.6 + 1e-8)  # asymptotes at beta
})

test_that("build_sigma_u returns a valid covariance matrix", {
  S <- build_sigma_u(c(2, 2.5, 3), rho = 0.4)
  expect_equal(dim(S), c(3L, 3L))
  expect_equal(diag(S), c(2, 2.5, 3))
  ev <- eigen(S, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that("build_sigma_u rejects bad variances and non-PD / out-of-range rho", {
  expect_error(build_sigma_u(c(2, 0, 3)), "positive")
  expect_error(build_sigma_u(c(2, 2.5, 3), rho = 1.2), "\\[-1, 1\\]")
  expect_error(build_sigma_u(c(2, 2.5, 3), rho = 1), "positive definite")
})

test_that("sigma_u_for_reliability validates var_x, ratios, rho", {
  expect_error(sigma_u_for_reliability(0.65, var_x = -5), "positive")
  expect_error(sigma_u_for_reliability(0.65, var_x = 5, ratios = c(2, 0, 3)),
               "positive")
  expect_error(sigma_u_for_reliability(0.65, var_x = 5, rho = 2), "\\[-1, 1\\]")
})

test_that("generate_aft_data rejects a negative x_sd", {
  expect_error(generate_aft_data(n = 50, n_val = 30, x_sd = -1, seed = 1),
               "non-negative")
})

test_that("generate_aft_data returns expected structure and reproduces with seed", {
  sim1 <- generate_aft_data(n = 200, n_val = 100, seed = 42)
  sim2 <- generate_aft_data(n = 200, n_val = 100, seed = 42)
  expect_identical(sim1$survival, sim2$survival)
  expect_identical(sim1$validation, sim2$validation)

  expect_equal(nrow(sim1$survival), 200L)
  expect_equal(nrow(sim1$validation), 100L)
  expect_true(all(c("T_obs", "delta", "X_true", "W1", "W2", "W3",
                    "V1", "V2", "V3", "V4") %in% names(sim1$survival)))
  expect_true(all(sim1$survival$delta %in% c(0L, 1L)))
})
