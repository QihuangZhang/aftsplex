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
