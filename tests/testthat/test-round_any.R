test_that("accuracy should be greater than 0", {
  expect_error(round_any(1.231, -0.005))
})

test_that("round works for positive x", {
  expect_true(round_any(1.231, 0.005) == 1.230)
})

test_that("round works for negative x", {
  expect_true(round_any(-1.231, 0.005) == -1.230)
})

test_that("floor works for positive x", {
  expect_true(round_any(1.231, 0.005, f = floor) == 1.230)
})

test_that("floor works for negative x", {
  expect_true(round_any(-1.231, 0.005, f = floor) == -1.235)
})

test_that("ceiling works for positive x", {
  expect_true(round_any(1.231, 0.005, f = ceiling) == 1.235)
})

test_that("ceiling works for negative x", {
  expect_true(round_any(-1.231, 0.005, f = ceiling) == -1.230)
})