test_that("extract_surrogate_cols finds W columns and errors when none match", {
  sim <- generate_aft_data(n = 50, n_val = 30, seed = 1)
  expect_equal(extract_surrogate_cols(sim$survival), c("W1", "W2", "W3"))
  expect_error(extract_surrogate_cols(data.frame(a = 1, b = 2)),
               "No surrogate columns matched")
})

test_that("fit_me_calibration returns alpha1 near 1 with large validation", {
  sim <- generate_aft_data(n = 50, n_val = 2000, seed = 1)
  cal <- fit_me_calibration(sim$validation)
  expect_length(cal$alpha0, 3L)
  expect_length(cal$alpha1, 3L)
  expect_equal(dim(cal$Sigma_e), c(3L, 3L))
  expect_true(all(abs(cal$alpha1 - 1) < 0.05))  # true alpha1 = 1 in DGP
})

test_that("gls_combine produces W_bar of correct length and positive sigma_w_sq", {
  sim <- generate_aft_data(n = 100, n_val = 500, seed = 1)
  cal <- fit_me_calibration(sim$validation)
  W <- as.matrix(sim$survival[, cal$W_cols])
  out <- gls_combine(W, cal)
  expect_length(out$W_bar, 100L)
  expect_gt(out$sigma_w_sq, 0)
  expect_equal(sum(out$omega), 1, tolerance = 1e-8)
})

test_that("fit_me_calibration guards degenerate validation inputs", {
  sim <- generate_aft_data(n = 50, n_val = 30, seed = 1)
  val_const <- sim$validation; val_const$X_true <- 5
  expect_error(fit_me_calibration(val_const), "constant or non-finite")
  expect_error(fit_me_calibration(sim$validation[1:2, ]), "at least 3")
  expect_error(fit_me_calibration(sim$validation, x_var = "nope"), "not found")
})

test_that("gls_combine errors informatively on degenerate calibration", {
  sim <- generate_aft_data(n = 80, n_val = 500, seed = 1)
  cal <- fit_me_calibration(sim$validation)
  W <- as.matrix(sim$survival[, cal$W_cols])

  W_bad <- W; W_bad[1, 1] <- NaN
  expect_error(gls_combine(W_bad, cal), "non-finite")

  cal0 <- cal; cal0$alpha1[2] <- 0
  expect_error(gls_combine(W, cal0), "no information")

  # collinear surrogates -> singular back-transformed covariance
  cal_sing <- cal; cal_sing$Sigma_e <- matrix(1, 3, 3)
  expect_error(gls_combine(W, cal_sing), "singular")
})
