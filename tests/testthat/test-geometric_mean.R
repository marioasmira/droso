test_that("geometric_mean doesn't work for negative values", {
  expect_error(geometric_mean(-50:0))
})

test_that("geometric_mean doesn't work for 0", {
  expect_error(geometric_mean(0))
})

test_that("geometric_mean works for positive values", {
  expect_true(geometric_mean(1:100) > 0)
})
