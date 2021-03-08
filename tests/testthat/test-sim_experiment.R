test_that("sim_experiment() validates inputs", {
  ref <- matrix(1L, nrow = 20, ncol = 10)
  ## reference
  expect_error(sim_experiment(1), "must be a matrix or sparseMatrix")
  expect_error(sim_experiment(1.5), "must be a matrix or sparseMatrix")
  expect_error(sim_experiment(T), "must be a matrix or sparseMatrix")
  expect_error(sim_experiment("A"), "must be a matrix or sparseMatrix")
  expect_error(sim_experiment(1:10), "must be a matrix or sparseMatrix")
  ## log2_fc
  expect_error(sim_experiment(ref, log2_fc = "A"), "must be a numeric")
  expect_error(sim_experiment(ref, log2_fc = T), "must be a numeric")
  expect_error(sim_experiment(ref, log2_fc = list(1, 2)), "must be a numeric")
  ## clones
  expect_error(sim_experiment(ref, n_clones = 1.5), "must be a integer")
  expect_error(sim_experiment(ref, n_clones = T), "must be a integer")
  expect_error(sim_experiment(ref, n_clones = "A"), "must be a integer")
  expect_error(sim_experiment(ref, n_clones = 1:10), "must be a integer scalar")
  expect_error(sim_experiment(ref, n_clones = -1), "must be a positive integer")
  expect_error(sim_experiment(ref, n_clones = 0), "must be a positive integer")
  ## genes
  expect_error(sim_experiment(ref, n_genes = 1.5), "must be a integer")
  expect_error(sim_experiment(ref, n_genes = T), "must be a integer")
  expect_error(sim_experiment(ref, n_genes = "A"), "must be a integer")
  expect_error(sim_experiment(ref, n_genes = 1:10), "must be a integer scalar")
  expect_error(sim_experiment(ref, n_genes = -1), "must be a positive integer")
  expect_error(sim_experiment(ref, n_genes = 0), "must be a positive integer")
})

test_that("sim_experiment() fails to sample more change than available", {
  ref <- matrix(1L, nrow = 20, ncol = 10)
  expect_error(
    sim_experiment(ref, log2_fc = rep(0, 20), n_genes = 10),
    "number of genes to be sampled for each clone exceed the total"
  )
})

test_that("sim_experiment() outputs the expected output", {
  nclones <- 10L
  ngenes <- 100L
  ncells <- 500L
  lfc <- rep(1:3, each = 5)
  ref <- matrix(1L, nrow = 500, ncol = ncells)
  res <- sim_experiment(ref, log2_fc = lfc, n_clones = nclones, n_genes = ngenes)
  expect_named(res, c("interactions", "clone", "count"))
  expect_equal(dim(res$count), c(ngenes, ncells))
  expect_equal(rowMeans(res$count), rep(1, ngenes), tolerance = 0.2)

  expect_equal(dim(res$clone), c(ncells, nclones))
  expect_equal(all(res$clone == 1 | res$clone == 0), TRUE)

  expect_named(res$interactions, c("clone", "gene", "log2fc"))
  expect_equal(nrow(res$interactions), nclones * 15)
  expect_equal(all(res$interactions$clone >= 1 & res$interactions$clone <= nclones), TRUE)
  expect_equal(all(res$interactions$gene >= 1 & res$interactions$clone <= ngenes), TRUE)
  expect_equal(res$interactions$log2fc, rep(lfc, nclones))
})
