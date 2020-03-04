
context("input_check")

test_that("input requirements are checked", {
  expect_error(input_check(data.frame(), data.frame()),
               regexp = "is not a matrix")
  expect_error(input_check(runif(5), data.frame()),
               regexp = "is not a matrix")
  expect_error(input_check(matrix(0, nrow = 5, ncol = 3),
                           matrix(0, nrow = 2, ncol = 2)),
               regexp = "is not a data frame")
  expect_error(input_check(matrix(0, nrow = 5, ncol = 3),
                           data.frame(matrix(0, nrow = 2, ncol = 3))),
               regexp = "not equal to")
  expect_error(input_check(matrix(0, nrow = 5, ncol = 3),
                           data.frame(matrix(0, nrow = 3, ncol = 3))),
               regexp = "does not have all of these name")
  expect_error(input_check(matrix(0, nrow = 5, ncol = 3),
                           data.frame(time = runif(3))),
               regexp = "does not have all of these name")
  expect_error(input_check(matrix(0, nrow = 5, ncol = 3),
                           data.frame(group = rep("A", 3))),
               regexp = "does not have all of these name")
  expect_error(input_check(cbind(1, 0, matrix(NA, nrow = 5)),
                           data.frame(group = rep("A", 3),
                                      time = runif(3))),
               regexp = "missing values")


})
