test_that("fit_aft_spline runs with and without covariates", {
  sim <- generate_aft_data(n = 300, n_val = 100, seed = 1)
  fit1 <- fit_aft_spline(sim$survival, x_var = "X_true",
                         covariates = c("V1","V2","V3","V4"))
  fit2 <- fit_aft_spline(sim$survival, x_var = "X_true",
                         covariates = character(0))
  expect_s3_class(fit1, "survreg")
  expect_s3_class(fit2, "survreg")
  # Adding covariates should change the coefficient count
  expect_gt(length(coef(fit1)), length(coef(fit2)))
})

test_that("predict_curve returns vector of correct length", {
  sim <- generate_aft_data(n = 300, n_val = 100, seed = 1)
  fit <- fit_aft_spline(sim$survival, x_var = "X_true",
                        covariates = c("V1","V2","V3","V4"))
  x_grid <- seq(7, 13, length.out = 40)
  lp <- predict_curve(fit, x_grid,
                      v_ref = c(V1=30, V2=30, V3=0, V4=0),
                      x_var = "X_true")
  expect_length(lp, 40L)
  expect_true(all(is.finite(lp)))
})

test_that("predict_curve accepts NULL v_ref when no covariates", {
  sim <- generate_aft_data(n = 300, n_val = 100, seed = 1)
  fit <- fit_aft_spline(sim$survival, x_var = "X_true",
                        covariates = character(0))
  x_grid <- seq(7, 13, length.out = 20)
  lp <- predict_curve(fit, x_grid, v_ref = NULL, x_var = "X_true")
  expect_length(lp, 20L)
})
