test_that("read_10x() input requirements", {
  expect_error(read_10x(1L), "Invalid `data_dir` argument")
  expect_error(read_10x(1.5), "Invalid `data_dir` argument")
  expect_error(read_10x(list(1, 2)), "Invalid `data_dir` argument")
  expect_error(read_10x(data.frame(x = 1:2)), "Invalid `data_dir` argument")
  expect_error(read_10x(c("A", "B")), "Invalid `data_dir` argument")
})

test_that("read_10x() loads matrix correctly", {
  data("sc_simulated")

  ## Create a temporary storage
  dir <- tempdir()
  mmfile <- file.path(dir, "matrix.mtx.gz")
  ftfile <- file.path(dir, "features.tsv.gz")
  bcfile <- file.path(dir, "barcodes.tsv.gz")
  Matrix::writeMM(sc_simulated, mmfile)
  writeLines(rownames(sc_simulated), ftfile)
  writeLines(colnames(sc_simulated), bcfile)

  mtx <- read_10x(dir)
  expect_equal(dim(mtx), dim(sc_simulated))
  expect_equal(rownames(mtx), rownames(sc_simulated))
  expect_equal(colnames(mtx), colnames(sc_simulated))
  expect_equal(sum(mtx - sc_simulated), 0)
})
