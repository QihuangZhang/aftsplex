test_that("compute_ise returns a non-negative ISE on a valid sorted grid", {
  x_grid <- seq(6, 14, length.out = 50)
  fits <- list(oracle = f_true(x_grid) - f_true(x_grid[1]),
               zero   = numeric(length(x_grid)))
  out <- compute_ise(x_grid, f_true, fits)
  expect_length(out$ise, 2L)
  expect_true(all(out$ise >= 0))
  expect_equal(unname(out$ise["oracle"]), 0)
})

test_that("compute_ise guards unsorted grids, short grids, and length mismatch", {
  x_grid <- seq(6, 14, length.out = 20)
  good <- f_true(x_grid) - f_true(x_grid[1])
  expect_error(compute_ise(rev(x_grid), f_true, list(a = rev(good))),
               "increasing order")
  expect_error(compute_ise(x_grid[1], f_true, list(a = good[1])),
               "at least 2 points")
  expect_error(compute_ise(x_grid, f_true, list(a = good[-1])),
               "length\\(x_grid\\)")
})
