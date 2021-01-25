test_that("col_div takes a matrix and apply a colwise division", {
  mat <- matrix(1L, nrow = 10, ncol = 5)
  res <- matrix(1/5, nrow = 10, ncol = 5)

  expect_error(col_div(mat, 1:2), "Incompatible dimensions")
  expect_equal(col_div(mat, rep(5, 5)), res)
})

test_that("col_div is compatible with sparse Matrix", {
  mat <- matrix(rep(0:4, each=10), ncol=5)
  sp1 <- as(mat, "dgCMatrix")
  sp2 <- as(mat, "dgTMatrix")
  res <- matrix(rep(c(0, 1/5, 2/5, 3/5, 4/5), each = 10), ncol=5)

  expect_equal(col_div(mat, rep(5, 5)), res)
  expect_equal(col_div(mat, rep(5, 5)), res)
  expect_equal(col_div(mat, rep(5, 5)), res)
})
