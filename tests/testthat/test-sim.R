ncell <- 100L
ngene <- 200L
dp <- sample(c(10436, 52619, 32060, 18266, 8254), ncell, replace = TRUE)
ps <- seq(0.01, 0.05, length.out = ngene)

check_sim <- function(sim) {
  sim_ps <- rowSums(sim) / sum(dp)
  expect_equal(is.matrix(sim), TRUE)
  expect_equal(ncol(sim), length(dp))
  expect_equal(nrow(sim), length(ps))
  expect_equal(rownames(sim), names(ps))
  expect_equal(colnames(sim), names(dp))
  expect_equal(mean((sim_ps - ps)**2), 0, tolerance = 1e-4)
}

test_that("sim_count() parameters type requirements", {
  unsupported_types <- list("a", FALSE, list(1), matrix(1, nrow = 2))

  test_type_requirement(unsupported_types, "must be a numeric",
    sim_count,
    gene_ps = ps, log2_fc = 0L
  )
  test_type_requirement(unsupported_types, "must be a numeric",
    sim_count,
    cell_dp = dp, log2_fc = 0L
  )
  test_type_requirement(unsupported_types, "must be a numeric",
    sim_count,
    cell_dp = dp, gene_ps = ps, log2_fc = 0L
  )
  test_type_requirement(unsupported_types[1:3], "must be a numeric",
    sim_count,
    cell_dp = dp, gene_ps = ps
  )
})

test_that("sim_count() parameters dimension requirements", {
  ## LOG2_FC SCALAR or N-Vector or MN-Matrix
  scalar_lfc <- 0L
  check_sim(sim_count(dp, ps, log2_fc = scalar_lfc))
  vector_lfc <- rep(0, ncell)
  check_sim(sim_count(dp, ps, log2_fc = vector_lfc))
  matrix_lfc <- matrix(0, nrow = ngene, ncol = ncell)
  check_sim(sim_count(dp, ps, log2_fc = matrix_lfc))
  ## Check size consistency
  expect_error(
    sim_count(dp, ps, log2_fc = rep(0, ncell - 5)),
    "Incompatible dimensions"
  )
  expect_error(
    sim_count(dp, ps, log2_fc = matrix(0, nrow = ngene, ncol = ncell - 5)),
    "Incompatible dimensions"
  )
  expect_error(
    sim_count(dp, ps, log2_fc = matrix(0, nrow = ngene - 5, ncol = ncell)),
    "Incompatible dimensions"
  )
})

test_that("sim_count() includes gene and cell names", {
  names(dp) <- sprintf("cell%03d", seq_along(dp))
  names(ps) <- sprintf("gene%03d", seq_along(ps))
  sim_mtx <- sim_count(dp, ps)
  expect_equal(colnames(sim_mtx), names(dp))
  expect_equal(rownames(sim_mtx), names(ps))
})

test_that("sim_count() produce desired log2 fold-change", {
  estimate_lfc <- function(sim, lfc) {
    lfc0 <- lfc == 0
    lfc1 <- lfc != 0
    ps0 <- rowSums(sim[, lfc0, drop = FALSE]) / sum(dp[lfc0])
    ps1 <- rowSums(sim[, lfc1, drop = FALSE]) / sum(dp[lfc1])
    ps1 / ps0
  }

  for (lfc in c(-2, -1, -.5, .5, 1, 2)) {
    lfc_vec <- rep(c(0, lfc), c(ncell - 10, 10))
    ## Vector lfc
    sim_mtx <- sim_count(dp, ps, log2_fc = lfc_vec)
    expect_equal(
      estimate_lfc(sim_mtx, lfc_vec), rep(2**lfc, ngene),
      tolerance = 0.2
    )
    ## Matrix lfc
    lfc_mtx <- t(replicate(ngene, sample(lfc_vec)))
    sim_mtx <- sim_count(dp, ps, log2_fc = lfc_mtx)
    est_lfc <- sapply(seq_len(ngene), function(i) {
      estimate_lfc(sim_mtx[i, , drop = FALSE], lfc_mtx[i, ])
    })
    expect_equal(est_lfc, rep(2**lfc, ngene), tolerance = 0.2)
  }
})

test_that("sim_count_from_data() requires matrix or sparseMatrix", {
  unsupported_types <- list(1, "a", FALSE, list(1), c(1, 2, 3))

  test_type_requirement(
    unsupported_types, "must be a matrix or sparseMatrix",
    sim_count_from_data
  )
})

test_that("sim_count_from_data() replicates similar data from base", {
  data("scSimulated")
  sim_mtx <- sim_count_from_data(scSimulated)
  expect_equal(dim(sim_mtx), dim(scSimulated))
  expect_equal(rowMeans(sim_mtx), rowMeans2(scSimulated), tolerance = 0.1)
  ## Small lfc change
  lfc_vec <- rep(c(0, 1), c(ncol(scSimulated) - 10, 10))
  sim_mtx <- sim_count_from_data(scSimulated, lfc_vec)
  expect_equal(dim(sim_mtx), dim(scSimulated))
  expect_equal(rowMeans(sim_mtx), rowMeans2(scSimulated), tolerance = 0.1)
})

test_that("sim_count_from_data() simulates n_genes as range", {
  sim_mtx <- sim_count_from_data(scSimulated, n_genes = 100L)
  expect_equal(ncol(sim_mtx), ncol(scSimulated))
  expect_equal(nrow(sim_mtx), 100L)

  exp_lmean <- log10(rowMeans2(scSimulated))
  sim_lmean <- log10(rowMeans(sim_mtx))
  expect_equal(min(10**sim_lmean), min(10**exp_lmean), tolerance = 0.05)
  expect_equal(max(10**sim_lmean), max(10**exp_lmean), tolerance = 0.05)
  exp_step <- diff(range(exp_lmean)) / 100L
  sim_step <- mean(diff(sim_lmean))
  expect_equal(sim_step, exp_step, tolerance = 0.05)
})
