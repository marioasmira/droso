# making sure the method without clipping the maximum value results in the same values
test_that("all values are equal", {
  x <- numeric()
  y <- numeric()
  for (i in 1:1000) {
    random_values <- c(
      18.5 + rnorm(1, 0, 5),
      1 + rnorm(1, 0, 1),
      34.5 + rnorm(1, 0, 5),
      -1.3 + rnorm(1, 0, 1),
      70 + rnorm(1, 0, 15)
    )
    f <- function(x) {
      squarelike(
        x,
        random_values[1],
        random_values[2],
        random_values[3],
        random_values[4],
        random_values[5]
      )
    }
    x <- c(x, unlist(integrate(f,-10, 40)[1]))

    y <- c(
      y,
      area_squarelike(
        -10,
        40,
        random_values[1],
        random_values[2],
        random_values[3],
        random_values[4],
        random_values[5]
      )
    )
    if (is.infinite(y[i]))
      print(random_values)
  }
  x <- unname(x)
  expect_equal(round(x), round(y))
})
