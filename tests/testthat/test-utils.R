unsupported <- list("a", FALSE, list(1), c(1, 2), list(1), matrix(1, nrow = 2))

test_that("seq_log() require numeric scalar `from` and `to`", {
  test_type_requirement(unsupported, "must be a numeric",
    seq_log, to = 10, length_out = 10L
  )
  test_type_requirement(unsupported, "must be a numeric",
    seq_log, from = 1, length_out = 10L
  )
})

test_that("seq_log() require integer scalar `length_out`", {
  test_type_requirement(c(unsupported, 10.5), "must be a integer scalar",
    seq_log, from = 1, to = 10
  )
  expect_error(seq_log(1, 10, -1), "must be a positive integer scalar")
  expect_error(seq_log(1, 10,  0), "must be a positive integer scalar")
})

test_that("seq_log() require numeric scalar `base`", {
  test_type_requirement(unsupported, "must be a numeric",
    seq_log, from = 1, to = 10, length_out = 10L
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
