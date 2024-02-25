test_that("correlation is greater or equal to 0", {
  expect_error(update_temperature(10, 100, -2, 25, 10, 1 / 20))
})

test_that("correlation is smaller or equal to 1", {
  expect_error(update_temperature(10, 100, 2, 25, 10, 1 / 20))
})