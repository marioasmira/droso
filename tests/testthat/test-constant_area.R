# making sure the method without clipping the maximum value results in the same values
test_that("all values are equal", {
  x <- numeric()
  y <- numeric()
  for(i in 1:1000){
    random_values <- c(
      18.5 + rnorm(1, 0, 10),
      1 + rnorm(1, 0, 10),
      34.5 + rnorm(1, 0, 10),
      -1.3 + rnorm(1, 0, 10),
      70 + rnorm(1, 0, 10)
    )
    f <- function(x) {
      ifelse(
        squarelike(
          x,
          random_values[1],
          random_values[2],
          random_values[3],
          random_values[4],
          random_values[5]
        ) < 0,
        0,
        squarelike(
          x,
          random_values[1],
          random_values[2],
          random_values[3],
          random_values[4],
          random_values[5]
        )
      )
    }
    x <- c(x, unlist(integrate(f, -10, 40)[1]))

    j <- function(x) {
      ifelse(
        squarelike(
          x,
          random_values[1],
          random_values[2],
          random_values[3],
          random_values[4],
          random_values[5]
        ) < 0,
        0,
        ifelse(
          squarelike(
            x,
            random_values[1],
            random_values[2],
            random_values[3],
            random_values[4],
            random_values[5]
          ) > random_values[5],
          random_values[5],
          squarelike(
            x,
            random_values[1],
            random_values[2],
            random_values[3],
            random_values[4],
            random_values[5]
          )
        )
      )
    }
    y <- c(y, unlist(integrate(j, -10, 40)[1]))
  }
  expect_equal(x, y)
})
