test_that("clone_nb input requirements", {
  X <- matrix(rbinom(50 * 1e3, 1, .1), ncol=50)
  Y <- matrix(rpois(1e3 * 1e3, 10), ncol=1e3)
  is <- sample(which(rowMeans(Y) > 1), 10)
  js <- sample(which(colSums(X) > 1), 10)

  expect_error(clone_nb(NULL, Y, X), "must be a data.frame")
  expect_error(clone_nb(integer(), Y, X), "must be a data.frame")
  expect_error(clone_nb(data.frame(gene=is, clone=js), integer(), X), "must be a matrix or sparseMatrix")
  expect_error(clone_nb(data.frame(gene=is, clone=js), Y, integer()), "must be a matrix or sparseMatrix")
  expect_error(clone_nb(data.frame(is, clone=js), y, X), "`pertubed` must contain")
  expect_error(clone_nb(data.frame(gene=is, js), y, X), "`pertubed` must contain")
})

test_that("clone_nb output is consistent with expected", {
  X <- matrix(rbinom(50 * 1e3, 1, .1), ncol=50)
  Y <- matrix(rpois(1e3 * 1e3, 10), ncol=1e3)
  d <- data.frame(
    gene=sample(which(rowMeans(Y) > 1), 10),
    clone=sample(which(colSums(X) > 1), 10),
    ignored = 0L)
  fit <- clone_nb(d, Y, X)
  output <- c(
    "n", "nonzero", "x1", "x0", "mu", "theta",
    "xb", "beta", "stderr", "z", "lfc", "pvalue", "qvalue")

  expect_type(fit, typeof(d))
  expect_named(fit, c(names(d), output), ignore.order = TRUE)
  expect_equal(mean(fit$pvalue < 0.05), 0.05, tolerance = 0.2)
  expect_equal(mean(fit$qvalue < 0.10), 0, tolerance = 0.2)
})

test_that("clone_nb add base computation when no test match criteria", {
  X <- matrix(rbinom(50 * 1e3, 1, .1), ncol=50)
  Y <- matrix(rpois(1e3 * 1e3, 10), ncol=1e3)
  d <- data.frame(
    gene=sample(which(rowMeans(Y) > 1), 10),
    clone=sample(which(colSums(X) > 1), 10),
    ignored = 0L)
  output <- c("n", "nonzero", "x1", "x0", "mu", "theta")
  ## forcing all tests to fail
  expect_warning(fit <- clone_nb(d, Y, X, min_n = 1e3), "passed the filter criteria")
  expect_type(fit, typeof(d))
  expect_named(fit, c(names(d), output), ignore.order = TRUE)
})

test_that("clone_nb filter criteria", {
  X <- matrix(rbinom(50 * 1e3, 1, .1), ncol=50)
  Y <- matrix(rpois(1e3 * 1e3, 10), ncol=1e3)
  X[, 1] <- sample(rep(0:1, c(1e3-1, 1)))  ## n < min_n [2]
  Y[5, ] <- rpois(1e3, .01)                     ## mu < min_tehta_umi [0.05]
  Y[7, ] <- rpois(1e3, .5)                      ## mu < min_x [1]
  fit <- clone_nb(
    data.frame(gene = c(5, 7, 10, 20), clone = c(25, 10, 7, 1)), Y, X)

  ## forcing all tests to fail
  expect_equal(is.na(fit$pvalue), c(T, T, F, T))
})
