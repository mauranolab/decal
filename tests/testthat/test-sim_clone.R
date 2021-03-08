test_that("sim_clone() and sim_clone_range() input requirements", {
  # size requirements
  expect_error(sim_clone(size = -1, n_cells = 10L), "must be a positive integer")
  expect_error(sim_clone(size = 0, n_cells = 10L), "must be a positive integer")
  expect_error(sim_clone(size = 1.5, n_cells = 10L), "must be a integer")
  expect_error(sim_clone(size = T, n_cells = 10L), "must be a integer")
  expect_error(sim_clone(size = "T", n_cells = 10L), "must be a integer")
  # n_cells requirements
  expect_error(sim_clone(size = 10L, n_cells = -1), "must be a positive integer")
  expect_error(sim_clone(size = 10L, n_cells = 0), "must be a positive integer")
  expect_error(sim_clone(size = 10L, n_cells = 1.5), "must be a integer")
  expect_error(sim_clone(size = 10L, n_cells = T), "must be a integer")
  expect_error(sim_clone(size = 10L, n_cells = "T"), "must be a integer")
  expect_error(sim_clone(size = 10L, n_cells = c(1, 2)), "must be a integer scalar")
  # min and max requirements
  expect_error(sim_clone_range(5, 10, min = -1), "must be a positive integer")
  expect_error(sim_clone_range(5, 10, min = 0), "must be a positive integer")
  expect_error(sim_clone_range(5, 10, min = 1.5), "must be a integer")
  expect_error(sim_clone_range(5, 10, min = T), "must be a integer")
  expect_error(sim_clone_range(5, 10, min = "T"), "must be a integer")
  expect_error(sim_clone_range(5, 10, min = c(1, 2)), "must be a integer scalar")

  expect_error(sim_clone_range(5, 10, max = -1), "must be a positive integer")
  expect_error(sim_clone_range(5, 10, max = 0), "must be a positive integer")
  expect_error(sim_clone_range(5, 10, max = 1.5), "must be a integer")
  expect_error(sim_clone_range(5, 10, max = T), "must be a integer")
  expect_error(sim_clone_range(5, 10, max = "T"), "must be a integer")
  expect_error(sim_clone_range(5, 10, max = c(1, 2)), "must be a integer scalar")
})

test_that("sim_clone() requires total size below n_cells", {
  expect_error(sim_clone(1:5, 10), "must be fewer or equal to `n_cells`")
})

test_that("sim_clone() produce the intended size without overlap", {
  mtx <- sim_clone(1:5, 20)
  expect_equal(nrow(mtx), 20)
  expect_equal(ncol(mtx), 5)
  expect_equal(colSums(mtx), 1:5)
  expect_equal(all(rowSums(mtx) <= 1), TRUE)
})

test_that("sim_clone_range() produce the intended size without overlap", {
  mtx <- sim_clone_range(5, 20, min = 2, max = 10)

  expect_equal(nrow(mtx), 20)
  expect_equal(ncol(mtx), 5)
  expect_equal(all(rowSums(mtx) <= 1), TRUE)
  expect_equal(all(colSums(mtx) <= 10 & colSums(mtx) >= 2), TRUE)
})

test_that("sim_clone_range() fails on unreachable clone size", {
  expect_error(
    sim_clone_range(5, 8, min = 2, max = 10),
    "Unable to sample clone size"
  )
})
