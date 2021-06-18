ncell <- 100L
ngene <- 200L
depth <- floor(rlnorm(ncell, 9, .5))
ratio <- seq(0.01, 0.05, length.out = ngene)
lfc <- rep(c(-2, -.5, 0, .5, 2), c(2, 2, 5, 2, 2))

estimate_lfc <- function(count, x) {
  r0 <- rowSums(count[, !x, drop = FALSE]) / sum(depth[!x])
  r1 <- rowSums(count[, x, drop = FALSE]) / sum(depth[x])
  return(log2(r1 / r0))
}

test_that("sim_experiment() outputs the expected output", {
  res <- sim_experiment(depth, ratio, lfc, nclones = 10L, min_n = 5, max_n = 20)
  expect_named(res, c("perturbations", "clone", "count"))
  ## Check count
  count <- res$count
  expect_equal(dim(count), c(ngene, ncell))
  expect_equal(var(rowSums(count) / sum(count) - ratio), 0, tolerance = 1e-4)
  ## Check clone
  clone <- res$clone
  sizes <- sapply(clone, length)
  expect_type(clone, "list")
  expect_equal(length(clone), 10L)
  expect_true(min(sizes) >= 5)
  expect_true(max(sizes) <= 20)
  expect_true(all(unlist(clone) %in% seq_len(ncell)))
  ## Check perturbations
  pert <- res$perturbations
  expect_named(pert, c("clone", "gene", "expected_lfc"))
  expect_equal(nrow(pert), 10 * length(lfc))
  expect_true(all(pert$clone >= 1 & pert$clone <= 10L))
  expect_true(all(pert$gene >= 1 & pert$gene <= ngene))
  ## Evaluate perturbation effect
  for (ci in seq_len(10L)) {
    xi <- seq_len(ncell) %in% clone[[ci]]
    gi <- pert$gene[pert$clone == ci]
    expect_equal(estimate_lfc(count[gi, ], xi), lfc, tolerance = 1)
  }
})

test_that("sim_experiment_from_data() outputs the expected output", {
  ref <- ratio %o% depth
  rat <- 10**(seq(log10(min(ratio)), log10(max(ratio)), length.out = 500L))

  res <- sim_experiment_from_data(ref, lfc, nclones = 10L, min_n = 5, max_n = 20, ngenes = 500L)
  expect_named(res, c("perturbations", "clone", "count"))
  ## Check count
  count <- res$count
  expect_equal(dim(count), c(500L, ncell))
  expect_equal(var(rowSums(count) / sum(count) - rat), 0, tolerance = 1e-3)
  ## Check clone
  clone <- res$clone
  sizes <- sapply(clone, length)
  expect_type(clone, "list")
  expect_equal(length(clone), 10L)
  expect_true(min(sizes) >= 5)
  expect_true(max(sizes) <= 20)
  expect_true(all(unlist(clone) %in% seq_len(ncell)))
  ## Check perturbations
  pert <- res$perturbations
  expect_named(pert, c("clone", "gene", "expected_lfc"))
  expect_equal(nrow(pert), 10 * length(lfc))
  expect_true(all(pert$clone >= 1 & pert$clone <= 10L))
  expect_true(all(pert$gene >= 1 & pert$gene <= 500L))
  ## Evaluate perturbation effect
  for (ci in seq_len(10L)) {
    xi <- seq_len(ncell) %in% clone[[ci]]
    gi <- pert$gene[pert$clone == ci]
    expect_equal(estimate_lfc(count[gi, ], xi), lfc, tolerance = 1)
  }
})

test_that("sim_experiment fails when request more perturbations than genes", {
  expect_error(
    sim_experiment(depth, seq(0.01, 0.05, length.out = 10), lfc = rep(0, 20)),
    "number of perturbations exceed the total number of genes"
  )
})

test_that("sim_experiment validates parametrs", {
  unsupported <- list("A", TRUE, list(1, 2))
  for (val in unsupported) {
    expect_error(sim_experiment(val, ratio), "must be a positive integer")
    expect_error(sim_experiment(depth, val), "must be numeric")
    expect_error(sim_experiment(depth, ratio, lfc = val), "must be numeric")
    expect_error(sim_experiment(depth, ratio, theta = val), "must be numeric")
    expect_error(sim_experiment(depth, ratio, nclones = val), "must be a positive integer scalar")
    expect_error(sim_experiment(depth, ratio, min_n = val), "must be a positive integer scalar")
    expect_error(sim_experiment(depth, ratio, max_n = val), "must be a positive integer scalar")
  }
  ## check positive integer
  for (val in c(-1, 0, 1.5)) {
    expect_error(sim_experiment(val, ratio), "must be a positive integer")
    expect_error(sim_experiment(depth, ratio, nclones = val), "must be a positive integer scalar")
    expect_error(sim_experiment(depth, ratio, min_n = val), "must be a positive integer scalar")
    expect_error(sim_experiment(depth, ratio, max_n = val), "must be a positive integer scalar")
  }
  ## check scalar
  expect_error(sim_experiment(depth, ratio, nclones = c(1, 2)), "must be a positive integer scalar")
  expect_error(sim_experiment(depth, ratio, min_n = c(1, 2)), "must be a positive integer scalar")
  expect_error(sim_experiment(depth, ratio, max_n = c(1, 2)), "must be a positive integer scalar")
  ## check rate
  expect_error(sim_experiment(depth, -1), "must be a ratio between 0 and 1")
  expect_error(sim_experiment(depth, 2), "must be a ratio between 0 and 1")
})


test_that("sim_experiment_from_data validates parametrs", {
  ref <- ratio %o% depth
  unsupported <- list("A", TRUE, list(1, 2))
  for (val in unsupported) {
    expect_error(sim_experiment_from_data(val), "must be a matrix")
    expect_error(sim_experiment_from_data(ref, lfc = val), "must be numeric")
    expect_error(sim_experiment_from_data(ref, theta = val), "must be numeric")
    expect_error(sim_experiment_from_data(ref, nclones = val), "must be a positive integer scalar")
    expect_error(sim_experiment_from_data(ref, min_n = val), "must be a positive integer scalar")
    expect_error(sim_experiment_from_data(ref, max_n = val), "must be a positive integer scalar")
    expect_error(sim_experiment_from_data(ref, ngenes = val), "must be a positive integer scalar")
  }
  ## check positive integer scalar
  for (val in list(-1, 0, 1.5, c(1, 2, 3))) {
    expect_error(sim_experiment_from_data(ref, nclones = val), "must be a positive integer scalar")
    expect_error(sim_experiment_from_data(ref, min_n = val), "must be a positive integer scalar")
    expect_error(sim_experiment_from_data(ref, max_n = val), "must be a positive integer scalar")
    expect_error(sim_experiment_from_data(ref, ngenes = val), "must be a positive integer scalar")
  }
  ## check matrix
  expect_error(sim_experiment_from_data(c(1, 2, 3, 4, 5)), "must be a matrix")
})
