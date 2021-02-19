test_that("as_index return x if x is integer and ref is not", {
  x <- sample(1:10)
  expect_equal(as_index(x), x)
  expect_equal(as_index(x, NULL), x)
  expect_equal(as_index(x, LETTERS[1:10]), x)
})

test_that("as_index returns index of x in ref if both are the same type", {
  i <- c(1, 2, 4, 6, 8, 10)

  expect_equal(as_index(LETTERS[i], LETTERS), i)
  expect_equal(as_index((10 + 1:10)[i], (10 + 1:10)), i)
})

test_that("as_index throw an error when x is not integer and ref is a different type", {
  expect_error(as_index(as_index(logical(1L), LETTERS)), "`x` must either be an integer or match `ref` type")
  expect_error(as_index(as_index(numeric(1L), LETTERS)), "`x` must either be an integer or match `ref` type")
  expect_error(as_index(as_index(character(1L), 1:10)), "`x` must either be an integer or match `ref` type")
})
