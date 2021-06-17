test_that("seq_log() parameters type requirements", {
  uns <- list(c(1, 2), list(1), FALSE, "C")
  for (i in seq_along(uns)) {
    expect_error(seq_log(uns[[i]], 10, 100), "must be a numeric scalar")
    expect_error(seq_log(1, uns[[i]], 100), "must be a numeric scalar")
    expect_error(seq_log(1, 10, 100, base=uns[[i]]), "must be a numeric scalar")
    expect_error(seq_log(1, 10, uns[[i]]), "must be a positive integer scalar")
  }
  expect_error(seq_log(1, 10, 0.5), "must be a positive integer scalar")
  expect_error(seq_log(1, 10, 0), "must be a positive integer scalar")
  expect_error(seq_log(1, 10, -1), "must be a positive integer scalar")

  expect_error(seq_log(10, 1, 10), "`from` must be smaller than `to`")
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
