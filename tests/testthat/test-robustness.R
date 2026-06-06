# Failure-mode handling: the package should fail loudly and informatively
# (clear error or warning) rather than crashing cryptically or returning
# silent NA. Motivated by a user hitting an all-NA bootstrap curve with no
# message and then an opaque error at plotting time.

test_that("collinear surrogates give a clear calibration error, not a Lapack message", {
  sim <- generate_aft_data(n = 120, n_val = 40, seed = 1)
  val <- sim$validation; sv <- sim$survival
  val$W2 <- val$W1; sv$W2 <- sv$W1   # perfectly collinear -> singular calibration
  expect_error(
    two_stage_bootstrap(
      survival = sv, validation = val, x_var = "X_true",
      covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30,V2=30,V3=0,V4=0),
      lambda = c(0.5, 1, 2), B = 3, R = 4, x_grid = seq(8, 12, length.out = 10)
    ),
    "calibration failed"
  )
})

test_that("a fully failing bootstrap warns and returns NA instead of crashing", {
  sim <- generate_aft_data(n = 120, n_val = 60, seed = 2)
  sv <- sim$survival; sv$delta <- 0L   # no events -> every fit fails
  w <- capture_warnings(
    out <- two_stage_bootstrap(
      survival = sv, validation = sim$validation, x_var = "X_true",
      covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30,V2=30,V3=0,V4=0),
      lambda = c(0.5, 1, 2), B = 3, R = 6, x_grid = seq(8, 12, length.out = 12)
    )
  )
  expect_match(w, "replicates failed", all = FALSE)
  expect_equal(out$R_effective, 0L)
  expect_true(all(is.na(out$lower)))
  expect_true(all(is.na(out$f_hat)))
})

test_that("simex returns graceful NA (not a '0 non-NA cases' error) when fits fail", {
  sim <- generate_aft_data(n = 120, n_val = 60, seed = 2)
  d <- sim$survival; d$delta <- 0L
  cal <- fit_me_calibration(sim$validation)
  d$W_bar <- gls_combine(as.matrix(d[, cal$W_cols]), cal)$W_bar
  w <- capture_warnings(
    out <- simex_aft_spline(
      d, x_var = "W_bar", sigma_w_sq = 1,
      covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30,V2=30,V3=0,V4=0),
      lambda = c(0.5, 1, 2), B = 3, x_grid = seq(8, 12, length.out = 10)
    )
  )
  expect_match(w, "extrapolation undefined", all = FALSE)
  expect_true(all(is.na(out$curve_simex)))   # graceful NA, no error
})

test_that("simex naive-fit failure (zero-row data) raises an informative error", {
  sim <- generate_aft_data(n = 120, n_val = 60, seed = 2)
  cal <- fit_me_calibration(sim$validation)
  d <- sim$survival[0, ]                       # no rows at all -> survreg cannot fit
  d$W_bar <- numeric(0)
  expect_error(
    suppressWarnings(simex_aft_spline(
      d, x_var = "W_bar", sigma_w_sq = 1,
      covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30,V2=30,V3=0,V4=0),
      lambda = c(0.5, 1, 2), B = 3, x_grid = seq(8, 12, length.out = 10)
    )),
    "naive AFT fit failed"
  )
})

test_that("plotting an all-NA curve stops with a clear message", {
  skip_if_not_installed("ggplot2")
  x <- seq(8, 12, length.out = 12)
  expect_error(
    suppressWarnings(
      plot_curves(x, truth = rep(0, 12), fits = list(SIMEX = rep(NA_real_, 12)))
    ),
    "Nothing to plot"
  )
})

test_that("plotting drops an all-NA series but keeps the valid ones", {
  skip_if_not_installed("ggplot2")
  x <- seq(8, 12, length.out = 12)
  expect_warning(
    p <- plot_curves(x, truth = rep(0, 12),
                     fits = list(SIMEX = seq(0, 1, length.out = 12),
                                 Naive = rep(NA_real_, 12))),
    "entirely NA"
  )
  expect_s3_class(p, "ggplot")
})

test_that("a clean bootstrap emits no failure warnings", {
  sim <- generate_aft_data(n = 250, n_val = 90, seed = 1)
  expect_no_warning(
    two_stage_bootstrap(
      survival = sim$survival, validation = sim$validation, x_var = "X_true",
      covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30,V2=30,V3=0,V4=0),
      lambda = c(0.5, 1, 2), B = 3, R = 6, x_grid = seq(7, 13, length.out = 20)
    )
  )
})
