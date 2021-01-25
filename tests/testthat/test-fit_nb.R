test_that("fit_nb return expected results under null", {
  Y <- matrix(rnbinom(1e5, 100, mu=10), ncol=1000)
  X <- matrix(rbinom(1000 * 50, 1, 0.1), ncol=50)
  is <- sample(which(rowMeans(Y) > 1), 1000, replace=TRUE)
  js <- sample(which(colSums(X) > 1), 1000, replace=TRUE)

  ft <- fit_nb(Y, X, 100, log(colSums(Y)), is, js)
  expect_equal(mean(ft$pvalue < 0.05), 0.05, tolerance = 0.1)
  expect_equal(mean(ft$qvalue < 0.10), 0, tolerance = 0.05)
})

test_that("fit_nb return expected under effect", {
  d <- rlnorm(1000, 10)
  x <- rbinom(1000, 1, 0.1)
  for (lfc in c(-5, -2, -1, 1, 2, 5)) {
    y <- rnbinom(1000, 100, mu = 5e-4 * d * 2^(lfc * x))
    f <- fit_nb(matrix(y, nrow=1), matrix(x, ncol=1), 100, log(d), 1, 1)
    expect_equal(f$lfc, lfc, tolerance = 0.25)
  }
})
