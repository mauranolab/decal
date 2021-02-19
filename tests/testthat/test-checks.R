test_that("validate_matrix()", {
  unsupported <- list(1, "a", T, list(1), data.frame(x = 1:2))
  test_type_requirement(
    unsupported, "must be a matrix or sparseMatrix",
    validate_matrix
  )
  matrix <- matrix(1L, nrow = 2, ncol = 2)
  expect_silent(validate_matrix(matrix))
  expect_silent(validate_matrix(as(matrix, "sparseMatrix")))
})

test_that("validate_adjacency_matrix()", {
  num_matrix <- matrix(c(-2, -1, 1, 2), nrow = 2, ncol = 2)
  lgl_matrix <- num_matrix < 0
  adj_matrix <- (num_matrix < 0) * 1L
  expect_error(
    validate_adjacency_matrix(num_matrix),
    "`num_matrix` must be a logical matrix or encoded as 0/1"
  )
  expect_silent(validate_adjacency_matrix(lgl_matrix))
  expect_silent(validate_adjacency_matrix(adj_matrix))
})

test_that("validate_dataframe()", {
  unsupported <- list(1, "a", T, list(1), matrix(1L, nrow = 2, ncol = 2))
  test_type_requirement(
    unsupported, "must be a data.frame like structure",
    validate_dataframe
  )
  expect_silent(validate_dataframe(data.frame(x = 1:10, y = 10:1)))
})

test_that("validate_numeric*()", {
  unsupported <- list("a", TRUE, list(1), matrix("a", nrow = 2, ncol = 2))
  num_scalar <- 1.5
  num_vector <- seq(1, 2, length.out = 10)
  num_matrix <- matrix(num_vector, nrow = 5)

  test_type_requirement(unsupported, "must be a numeric", validate_numeric)
  expect_silent(validate_numeric(num_scalar))
  expect_silent(validate_numeric(num_vector))
  expect_silent(validate_numeric(num_matrix))

  expect_error(validate_numeric_vector(num_matrix), "must be a numeric vector")
  expect_silent(validate_numeric_vector(num_scalar))
  expect_silent(validate_numeric_vector(num_vector))

  expect_error(validate_numeric_scalar(num_matrix), "must be a numeric scalar")
  expect_error(validate_numeric_scalar(num_vector), "must be a numeric scalar")
  expect_silent(validate_numeric_scalar(num_scalar))
})

test_that("validate*_integer_scalar()", {
  unsupported <- list("a", TRUE, list(1), 1:5, matrix(1L, nrow = 2, ncol = 2))
  pos_scalar <- 1L
  neg_scalar <- -1L
  zero_scalar <- 0

  test_type_requirement(
    unsupported, "must be a integer scalar",
    validate_integer_scalar
  )
  expect_silent(validate_integer_scalar(pos_scalar))
  expect_silent(validate_integer_scalar(neg_scalar))
  expect_silent(validate_integer_scalar(zero_scalar))

  expect_silent(validate_positive_integer_scalar(pos_scalar))
  expect_error(
    validate_positive_integer_scalar(neg_scalar),
    "must be a positive integer scalar"
  )
  expect_error(
    validate_positive_integer_scalar(zero_scalar),
    "must be a positive integer scalar"
  )
})
