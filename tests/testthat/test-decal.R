test_that("decal recognize named indexes", {
  cname <- sprintf("cell-%03d", seq_len(100))
  rname <- sprintf("gene-%03d", seq_len(100))
  gname <- sprintf("clone-%03d", seq_len(5))
  cols <- c(
    "gene", "clone", "n0", "n1", "x0", "x1", "mu", "theta",
    "xb", "z", "lfc", "pvalue", "p_adjusted"
  )

  count <- matrix(rpois(100 * 100, 10), ncol=100, dimnames = list(rname, cname))
  clone <- split(seq_len(100), gname[sample(1:5, 100, replace = TRUE)])
  clone <- sapply(clone, function(x) cname[x])
  perturbations <- data.frame(
    gene = rname[sample(which(rowMeans(count) > 1), 20, replace=TRUE)],
    clone = gname[sample(1:5, 20, replace=TRUE)]
  )

  result <- decal(perturbations, count, clone)
  expect_type(result, typeof(perturbations))
  expect_equal(nrow(result), nrow(perturbations))
  expect_named(result, cols, ignore.order = TRUE)
  expect_equal(result$mu, rowMeans(count)[perturbations$gene],
               ignore_attr = TRUE)
})

test_that("decal validates perturbation format", {
  count <- matrix(rpois(100 * 100, 10), ncol=100)
  clone <- split(seq_len(100), sample(1:5, 100, replace = TRUE))
  perturbations <- data.frame(
    x=sample(which(rowMeans(count) > 1), 5, replace=TRUE),
    y=sample(1:5, 5, replace=TRUE)
  )
  expect_error(decal(integer(), count, clone), "must be a data.frame")
  expect_error(decal(matrix(), count, clone), "must be a data.frame")
  expect_error(decal(perturbations, count, clone),
               "must contain all required columns")
})

test_that("decal validates count as matrix", {
  count <- matrix(rpois(100 * 100, 10), ncol=100)
  clone <- split(seq_len(100), sample(1:5, 100, replace = TRUE))
  perturbations <- data.frame(
    gene=sample(which(rowMeans(count) > 1), 5, replace=TRUE),
    clone=sample(1:5, 5, replace=TRUE)
  )
  expect_error(decal(perturbations, integer(), clone),
               "must be a matrix")
  expect_error(decal(perturbations, data.frame(), clone),
               "must be a matrix")
})

test_that("decal validates clone as list", {
  count <- matrix(rpois(100 * 100, 10), ncol=100)
  clone <- split(seq_len(100), sample(1:5, 100, replace = TRUE))
  perturbations <- data.frame(
    gene=sample(which(rowMeans(count) > 1), 5, replace=TRUE),
    clone=sample(1:5, 5, replace=TRUE)
  )
  expect_error(decal(perturbations, count, integer()),
               "must be a list")
  expect_error(decal(perturbations, count, matrix()),
               "must be a list")
})

test_that("decal output is consistent with perturbations", {
  count <- matrix(rpois(100 * 100, 10), ncol=100)
  clone <- split(seq_len(100), sample(1:5, 100, replace = TRUE))
  perturbations <- data.frame(
    x=sample(which(rowMeans(count) > 1), 100, replace=TRUE),
    y=sample(1:5, 10, replace=TRUE)
  )
  expected_cols <- c("x", "y", "n0", "n1", "x0", "x1", "mu", "theta",
                     "xb", "z", "lfc", "pvalue", "p_adjusted")

  result <- decal(perturbations, count, clone, gene_col="x", clone_col="y")
  expect_type(result, typeof(perturbations))
  expect_equal(nrow(result), nrow(perturbations))
  expect_named(result, expected_cols, ignore.order = TRUE)
  expect_equal(mean(result$pvalue < 0.05), 0.05, tolerance = 0.1)
  expect_equal(mean(result$p_adjusted < 0.10), 0, tolerance = 0.1)
})

test_that("decal skips tests when don't match criterias", {
  X <- split(seq_len(100), sample(1:50, 100, replace = TRUE))
  Y <- matrix(rpois(1e3 * 1e3, 3), ncol = 1e3)
  d <- data.frame(
    gene = sample(which(rowMeans(Y) > 1), 10),
    clone = sample(which(sapply(X, length) > 1), 10)
  )
  expected_cols <- c("n0", "n1", "x0", "x1", "mu", "theta",
                     "xb", "z", "lfc", "pvalue", "p_adjusted")
  ## forcing tests to fail
  expect_warning(
    fit <- decal(d, Y, X, min_n = 1E3),
    "No gene & clone perturbation matched the requirements")
  expect_type(fit, typeof(d))
  expect_named(fit, c(names(d), expected_cols), ignore.order = TRUE)
  expect_equal(sum(is.na(fit$pvalue)), 10)
})
