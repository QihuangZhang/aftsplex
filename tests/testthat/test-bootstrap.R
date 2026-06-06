test_that("two_stage_bootstrap returns correct shape with tiny R, B", {
  sim <- generate_aft_data(n = 200, n_val = 80, seed = 1)
  x_grid <- seq(7, 13, length.out = 20)
  set.seed(1)
  boot <- two_stage_bootstrap(
    survival   = sim$survival,
    validation = sim$validation,
    x_var      = "X_true",
    covariates = c("V1","V2","V3","V4"),
    v_ref      = c(V1 = 30, V2 = 30, V3 = 0, V4 = 0),
    df         = 4,
    lambda     = c(0.5, 1, 2), B = 3, R = 4,
    x_grid     = x_grid
  )
  expect_length(boot$f_hat, 20L)
  expect_equal(dim(boot$curves), c(20L, 4L))
  expect_equal(length(boot$lower), 20L)
  expect_equal(length(boot$upper), 20L)
  expect_equal(length(boot$lower_basic), 20L)
  expect_equal(length(boot$f_bc), 20L)
  expect_equal(boot$R_effective, 4L)
  # The two CI definitions reflect around f_hat
  expect_equal(boot$lower_basic, 2 * boot$f_hat - boot$upper)
  expect_equal(boot$upper_basic, 2 * boot$f_hat - boot$lower)
  expect_equal(boot$f_bc,        2 * boot$f_hat - boot$bs_mean)
})

test_that("nested = FALSE holds the Phase-1 calibration fixed", {
  sim <- generate_aft_data(n = 200, n_val = 80, seed = 1)
  args <- list(
    survival = sim$survival, validation = sim$validation, x_var = "X_true",
    covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30, V2=30, V3=0, V4=0),
    lambda = c(0.5, 1), B = 3, R = 6, x_grid = seq(7, 13, length.out = 15)
  )
  set.seed(1); single <- do.call(two_stage_bootstrap, c(args, nested = FALSE))
  set.seed(1); double <- do.call(two_stage_bootstrap, c(args, nested = TRUE))
  # sigma_w_sq / omega depend only on the calibration, so a fixed calibration
  # gives identical values across replicates; the double bootstrap varies them.
  expect_equal(stats::sd(single$sigma_w_sq_boot), 0, tolerance = 1e-12)
  expect_gt(stats::sd(double$sigma_w_sq_boot), 0)
  expect_equal(dim(single$curves), dim(double$curves))
})

test_that("x_ref re-anchors the curve without changing its shape", {
  sim <- generate_aft_data(n = 200, n_val = 80, seed = 1)
  x_grid <- seq(7, 13, length.out = 21)
  args <- list(
    survival = sim$survival, validation = sim$validation, x_var = "X_true",
    covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30, V2=30, V3=0, V4=0),
    lambda = c(0.5, 1), B = 3, R = 4, x_grid = x_grid
  )
  set.seed(2); base <- do.call(two_stage_bootstrap, args)
  set.seed(2); ref  <- do.call(two_stage_bootstrap, c(args, x_ref = x_grid[11]))
  # curve passes through zero at x_ref ...
  expect_equal(unname(ref$f_hat[11]), 0, tolerance = 1e-10)
  # ... and the shape (first differences) is unchanged.
  expect_equal(diff(ref$f_hat), diff(base$f_hat), tolerance = 1e-10)
})

test_that("support_probs flags out-of-support grid points", {
  sim <- generate_aft_data(n = 200, n_val = 80, seed = 1)
  x_grid <- seq(2, 18, length.out = 25)   # deliberately exceeds the support
  set.seed(3)
  expect_warning(
    boot <- two_stage_bootstrap(
      survival = sim$survival, validation = sim$validation, x_var = "X_true",
      covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30, V2=30, V3=0, V4=0),
      lambda = c(0.5, 1), B = 3, R = 4, x_grid = x_grid,
      support_probs = c(0.05, 0.95)
    ),
    "extrapolations"
  )
  expect_length(boot$in_support, length(x_grid))
  expect_type(boot$in_support, "logical")
  expect_false(boot$in_support[1])             # low end out of support
  expect_true(any(boot$in_support))            # middle in support
})

test_that("workers > 1 is reproducible under a fixed seed", {
  skip_if_not_installed("future.apply")
  sim <- generate_aft_data(n = 200, n_val = 80, seed = 1)
  args <- list(
    survival = sim$survival, validation = sim$validation, x_var = "X_true",
    covariates = c("V1","V2","V3","V4"), v_ref = c(V1=30, V2=30, V3=0, V4=0),
    lambda = c(0.5, 1), B = 3, R = 4, x_grid = seq(7, 13, length.out = 12),
    workers = 2
  )
  set.seed(7); a <- do.call(two_stage_bootstrap, args)
  set.seed(7); b <- do.call(two_stage_bootstrap, args)
  expect_equal(a$curves, b$curves)
})
