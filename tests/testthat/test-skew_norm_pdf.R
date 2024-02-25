# normal_pdf
test_that("normal pdf returns only positive values", {
  expect_true(all(normal_pdf(rnorm(n = 1000, mean = 0, sd = 1)) > 0))
})

test_that("normal pdf returns error if no argument", {
  expect_error(normal_pdf())
})

# normal_cdf
test_that("normal cdf returns only positive values", {
  expect_true(all(normal_cdf(rnorm(n = 1000, mean = 0, sd = 1)) > 0))
})

test_that("normal cdf returns error if no argument", {
  expect_error(normal_cdf())
})

# skew_normal_pdf

test_that("skew_normal_pdf requires positive sd", {
  expect_error(skew_normal_pdf(5, 10, -2, 3))
})

test_that("skew_normal_pdf requires all arguments present", {
  expect_error(skew_normal_pdf(, 10, -2, 3))
})

test_that("skew_normal_pdf requires all arguments present", {
  expect_error(skew_normal_pdf(numeric(), 10, -2, 3))
})