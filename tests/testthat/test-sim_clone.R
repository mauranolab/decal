test_that("sim_clone() and sim_clone_range() input requirements", {
  vars <- list(-1, 0, 0.5, TRUE, "C")
  # clone_size requirements
  for (i in seq_along(vars)) {
    expect_error(sim_clone(vars[[i]], 10L), "must be a positive integer")
    expect_error(sim_clone(10L, vars[[i]]), "must be a positive integer")
    expect_error(sim_clone_range(vars[[i]], 10L), "must be a positive integer")
    expect_error(sim_clone_range(2, 10L, min = vars[[i]]), "must be a positive integer")
    expect_error(sim_clone_range(2, 10L, max = vars[[i]]), "must be a positive integer")
  }
  expect_error(sim_clone_range(2, 10L, min = 10, max = 2),
               "Minimum number of cells per clone must be less than the maximum")
  expect_error(sim_clone_range(10, 10L),
               "Unable to generate clones with current `ncells` limit")
})

test_that("sim_clone() produces expected clone", {
  clone <- sim_clone(2:10, 100L)
  expect_type(clone, "list")
  expect_equal(length(clone), 9)
  expect_equal(sapply(clone, length), 2:10)
  expect_true(all(unlist(clone) %in% seq_len(100)))
})

test_that("sim_clone_range() produces expected clones", {
  clone <- sim_clone_range(10, 100L, 2, 10)
  clone_size <- sapply(clone, length)
  expect_type(clone, "list")
  expect_equal(length(clone), 10)
  expect_true(min(clone_size) >= 2)
  expect_true(max(clone_size) <= 10)
  expect_true(all(unlist(clone) %in% seq_len(100)))
})
