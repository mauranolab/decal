ncell <- 100L
ngene <- 200L
depth <- floor(rlnorm(ncell, 9.5, .5))
ratio <- seq(0.01, 0.05, length.out = ngene)

test_that("sim_count produce the expected output", {
  res <- sim_count(depth, ratio)
  expect_true(is.matrix(res))
  expect_equal(ncol(res), ncell)
  expect_equal(nrow(res), ngene)
  expect_equal(rowSums(res) / sum(res), ratio, tolerance = 0.05)
})

test_that("sim_count_from_data produce the expected output", {
  mu <- ratio %o% depth
  ref <- matrix(rnbinom(ncell * ngene, size = 100L, mu = mu),
                ncol = ncell, nrow = ngene)
  res <- sim_count_from_data(ref)
  expect_true(is.matrix(res))
  expect_equal(ncol(res), ncell)
  expect_equal(nrow(res), ngene)
  expect_equal(var(rowSums(res) / sum(res) - ratio), 0, tolerance = 1e-4)
})

test_that("sim_count_from_data with ngenes produce the expected output", {
  ## NOTE: THAT FOR SMALL N THE INDEPENDENCE DOESN'T HOLD AND RATIO DEVIATES
  ## FROM EXPECTED
  n   <- 50
  ref <- ratio %o% depth
  rat <- 10 ** (seq(log(min(ratio), 10), log(max(ratio), 10), length.out = n))
  res <- sim_count_from_data(ref, ngenes = n)

  expect_true(is.matrix(res))
  expect_equal(ncol(res), ncell)
  expect_equal(nrow(res), n)
  expect_equal(var(rowSums(res) / sum(res) - rat), 0, tolerance = 1e-4)
})

test_that("sim_clone checks parameters", {
  unsupported <- list("C", TRUE, list(1, 2))
  for (val in unsupported) {
    expect_error(sim_count(val, ratio), "must be a positive integer")
    expect_error(sim_count(depth, val), "must be numeric")
    expect_error(sim_count(depth, ratio, theta = val), "must be numeric")
    expect_error(sim_count(depth, ratio, lfc = val), "must be numeric")
  }
  ## check is positive integer
  expect_error(sim_count(depth = -1, ratio), "must be a positive integer")
  expect_error(sim_count(depth = 0, ratio), "must be a positive integer")
  expect_error(sim_count(depth = 1.2, ratio), "must be a positive integer")
  ## check is a rate
  expect_error(sim_count(depth, -1), "must be a ratio between 0 and 1")
  expect_error(sim_count(depth, 1.2), "must be a ratio between 0 and 1")
})

test_that("sim_clone_from_data checks parameters", {
  ref <- ratio %o% depth
  unsupported <- list("C", TRUE, list(1, 2))
  for (val in unsupported) {
    expect_error(sim_count_from_data(val), "must be a matrix")
    expect_error(sim_count_from_data(ref, lfc=val), "must be numeric")
    expect_error(sim_count_from_data(ref, theta=val), "must be numeric")
    expect_error(sim_count_from_data(ref, ngenes=val), "must be a positive integer")
  }
  ## check only takes matrices
  expect_error(sim_count_from_data(1), "must be a matrix")
  expect_error(sim_count_from_data(c(1, 2)), "must be a matrix")
  ## check is positive integer
  expect_error(sim_count_from_data(ref, ngenes=-1), "must be a positive integer")
  expect_error(sim_count_from_data(ref, ngenes=0), "must be a positive integer")
  expect_error(sim_count_from_data(ref, ngenes=1.2), "must be a positive integer")
})

test_that("sim_count check dimensions requirements", {
  expect_error(
    sim_count(depth, ratio, theta = rep(100, 5)),
    "`theta` must be a scalar or vector of same length as `ratio`"
  )
  expect_error(
    sim_count(depth, ratio, lfc = rep(0, 5)),
    "`lfc` must be a scalar, vector of same length as `depth` or a matrix"
  )
  ## Check lfc matrix dimensions
  expect_error(sim_count(depth, ratio, lfc = matrix(0, nrow=10, ncol=10)),
               "Incompatible dimensions")
  expect_error(sim_count(depth, ratio, lfc = matrix(0, nrow=ngene, ncol=10)),
               "Incompatible dimensions")
  expect_error(sim_count(depth, ratio, lfc = matrix(0, nrow=10, ncol=ncell)),
               "Incompatible dimensions")
})

test_that("sim_count includes gene and cell names", {
  names(depth) <- sprintf("cell%03d", seq_along(depth))
  names(ratio) <- sprintf("gene%03d", seq_along(ratio))
  res <- sim_count(depth, ratio)
  expect_equal(colnames(res), names(depth))
  expect_equal(rownames(res), names(ratio))
})

test_that("sim_count_from_data preserves gene and cell names", {
  names(depth) <- sprintf("cell%03d", seq_along(depth))
  names(ratio) <- sprintf("gene%03d", seq_along(ratio))
  ref <- ratio %o% depth
  res <- sim_count_from_data(ref)
  expect_equal(colnames(res), names(depth))
  expect_equal(rownames(res), names(ratio))
})

test_that("sim_count produces the desired log fold-change", {
  estimate_lfc <- function(mat, x) {
    r0 <- rowSums(mat[, x == 0, drop=FALSE]) / sum(depth[x == 0])
    r1 <- rowSums(mat[, x == 1, drop=FALSE]) / sum(depth[x == 1])
    log2(r1 / r0)
  }

  x <- sample(rep(0:1, c(ncell - 20, 20)))
  for (lfc in c(-2, -.5, 0, .5, 2)) {
    ## Vectorized LFC
    lfc_vec <- x * lfc
    res <- sim_count(depth, ratio, lfc_vec)
    expect_equal(estimate_lfc(res, x), rep(lfc, ngene), tolerance = 0.1)
    ## Matrix LFC
    lfc_mat <- replicate(ngene, lfc_vec)
    res <- sim_count(depth, ratio, lfc_vec)
    expect_equal(estimate_lfc(res, x), rep(lfc, ngene), tolerance = 0.1)
  }
})
