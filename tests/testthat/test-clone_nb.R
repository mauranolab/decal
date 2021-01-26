test_that("clone_nb input validation", {
  X <- matrix(rbinom(50 * 1e2, 1, .1), ncol=50)
  Y <- matrix(rpois(1e2 * 1e2, 10), ncol=1e2)
  is <- sample(which(rowMeans(Y) > 1), 10)
  js <- sample(which(colSums(X) > 1), 10)

  expect_error(clone_nb(integer(), Y, X), "must be a data.frame")
  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, gene_col=logical()),
    "must be either a column name or index")
  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, clone_col=logical()),
    "must be either a column name or index")
  expect_error(clone_nb(data.frame(is, clone=js), Y, X), "is not a column name")
  expect_error(clone_nb(data.frame(gene=is, js), Y, X), "is not a column name")

  expect_error(clone_nb(data.frame(gene=is, clone=js), integer(), X),
    "must be a matrix or sparseMatrix")
  expect_error(clone_nb(data.frame(gene=is, clone=js), Y, integer()),
    "must be a matrix or sparseMatrix")
  expect_error(clone_nb(data.frame(gene=is, clone=js), Y, Y),
    "must be a logical matrix or encoded as 0/1")
  expect_error(clone_nb(data.frame(gene=is, clone=js), Y, X[1:10,]),
    "have incompatible dimensions")

  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, min_x = logical()),
    "must be a scalar numeric value")
  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, min_n = logical()),
    "must be a scalar numeric value")
  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, theta_min_mu = logical()),
    "must be a scalar numeric value")
  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, theta_n = logical()),
    "must be a scalar numeric value")

  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, min_x = integer(2L)),
    "must be a scalar numeric value")
  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, min_n = integer(2L)),
    "must be a scalar numeric value")
  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, theta_min_mu = integer(2L)),
    "must be a scalar numeric value")
  expect_error(
    clone_nb(data.frame(gene=is, clone=js), Y, X, theta_n = integer(2L)),
    "must be a scalar numeric value")
})

test_that("clone_nb result is consistent with input variation", {
  X <- matrix(rbinom(50 * 1e2, 1, .1), ncol=50)
  Y <- matrix(rpois(1e2 * 1e2, 10), ncol=1e2)
  dt <- data.frame(
    gene = sample(which(rowMeans(Y) > 1), 10L),
    clone = sample(which(colSums(X) > 1), 10L)
  )
  ft <- clone_nb(dt, Y, X)

  expect_equal(clone_nb(dt, Y, X == 1), ft)
  expect_equal(clone_nb(dt, Y, X, gene_col = 1L), ft)
  expect_equal(clone_nb(dt, Y, X, clone_col = 2L), ft)
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
