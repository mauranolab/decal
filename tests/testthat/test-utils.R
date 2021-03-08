unsupported <- list("a", FALSE, list(1), c(1, 2), list(1), matrix(1, nrow = 2))

test_that("seq_log requires numeric scalar `from` and `to`", {
  test_type_requirement(unsupported, "must be a numeric",
    seq_log,
    to = 10, length_out = 10L
  )
  test_type_requirement(unsupported, "must be a numeric",
    seq_log,
    from = 1, length_out = 10L
  )
})

test_that("seq_log() require integer scalar `length_out`", {
  test_type_requirement(c(unsupported, 10.5), "must be a integer",
    seq_log,
    from = 1, to = 10
  )
  expect_error(seq_log(1, 10, -1), "must be a positive integer scalar")
  expect_error(seq_log(1, 10, 0), "must be a positive integer scalar")
})

test_that("seq_log() require numeric scalar `base`", {
  test_type_requirement(unsupported, "must be a numeric",
    seq_log,
    from = 1, to = 10, length_out = 10L
  )
})

test_that("seq_log() produces an equally spaced sequence in log scale", {
  expect_equal(seq_log(1, 2**2, length_out = 3, base = 2), 2**(0:2))
  expect_equal(seq_log(1, 10**2, length_out = 3, base = 10), 10**(0:2))
  expect_equal(seq_log(1, exp(2), length_out = 3, base = exp(1)), exp(0:2))
})

test_that("seq_log() produces desired length", {
  expect_equal(length(seq_log(1, 10, length_out = 2)), 2)
  expect_equal(length(seq_log(1, 10, length_out = 10)), 10)
  expect_equal(length(seq_log(1, 10, length_out = 100)), 100)
  len_one_seq <- seq_log(1, 10, length_out = 1)
  expect_equal(length(len_one_seq), 1)
})

test_that("seq_log() min and max matches from and to", {
  len_one_seq <- seq_log(1, 10, length_out = 1)
  expect_equal(len_one_seq, 1)
  len_ten_seq <- seq_log(1, 10, length_out = 10)
  expect_equal(min(len_ten_seq), 1)
  expect_equal(max(len_ten_seq), 10)
})

## as_clone_matrix
test_that("as_clone_matrix() requires same length index vector", {
  ## clones type-check
  test_type_requirement(list(1.5, FALSE), "must be a integer or character",
    as_clone_matrix,
    cells = 1
  )
  expect_error(as_clone_matrix(matrix(1, nrow = 2), cells = 1:2), "must be a integer or character")
  ## cells type check
  test_type_requirement(list(1.5, FALSE), "must be a integer or character",
    as_clone_matrix,
    clones = 1
  )
  expect_error(as_clone_matrix(matrix(1, nrow = 2), clones = 1:2), "must be a integer or character")
  ## length check
  expect_error(as_clone_matrix(1, 1:2), "must have the same length")
})

# test_that("as_clone_matrix() accepts both vector or list `cells`", {
#   clones_vector <- rep(1:10, 3)
#   cells_vector <- 1:30
#   mtx_vector <- as_clone_matrix(clones_vector, cells_vector)
#   mtx_list <- as_clone_matrix(1:10, split(cells_vector, clones_vector))
#   expect_equal(mtx_vector, mtx_list)
# })
#
# test_that("as_clone_matrix() fills blank with integer indexes", {
#   mtx <- as_clone_matrix(c(1, 1, 1, 3, 3, 4), c(1, 2, 3, 6, 7, 10))
#   expect_equal(nrow(mtx), 10)
#   expect_equal(ncol(mtx), 4)
#   expect_equal(sum(mtx[, 2]), 0)
#   expect_equal(sum(mtx[4, ]), 0)
#   expect_equal(sum(mtx[9, ]), 0)
# })
#
# test_that("as_clone_matrix() orders character names", {
#   cnm <- letters[c(1, 3, 4)]
#   rnm <- LETTERS[c(1, 2, 3, 6, 7, 10)]
#   mtx <- as_clone_matrix(cnm[c(3, 2, 2, 1, 1, 1)], rnm[c(6, 5, 4, 3, 2, 1)])
#   expect_equal(nrow(mtx), length(rnm))
#   expect_equal(ncol(mtx), length(cnm))
#   expect_equal(rownames(mtx), rnm)
#   expect_equal(colnames(mtx), cnm)
# })
