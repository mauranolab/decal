decal_cols <- c("n0", "n1", "x0", "x1", "mu", "theta", "xb", "z", "lfc", "pvalue", "p_adjusted")

build_dat <- function(ngene = 100, ncell = 100, nclone = 5, ntest = 20, named = FALSE) {
  count <- matrix(rpois(ngene * ncell, 10), ncol = ncell)
  clone <- split(
    seq_len(ncell),
    sample(seq_len(nclone), ncell, replace = TRUE)
  )
  perturbations <- data.frame(
    gene = sample(which(rowMeans(count) > 1), ntest, replace = TRUE),
    clone = sample(1:5, ntest, replace = TRUE)
  )
  if (named == TRUE) {
    cname <- sprintf("cell-%03d", seq_len(100))
    rname <- sprintf("gene-%03d", seq_len(100))
    gname <- sprintf("clone-%03d", seq_len(5))

    dimnames(count) <- list(rname, cname)
    clone <- sapply(clone, function(x) cname[x])
    names(clone) <- gname
    perturbations$gene <- rname[perturbations$gene]
    perturbations$clone <- gname[perturbations$clone]
  }
  return(list(pert = perturbations, count = count, clone = clone))
}

check_decal <- function(dat, result) {
  ## Check result consistency
  cols <- c(colnames(dat$pert), decal_cols)
  gene_mean <- rowMeans(dat$count)
  clone_size <- sapply(dat$clone, length)
  ## Check values
  expect_type(result, typeof(dat$pert))
  expect_equal(nrow(result), nrow(dat$pert))
  expect_named(result, cols, ignore.order = TRUE)
  expect_equal(result$mu, gene_mean[dat$pert$gene], ignore_attr = TRUE)
  expect_equal(result$n1, clone_size[dat$pert$clone], ignore_attr = TRUE)
}

test_that("decal produces the expected result", {
  dat <- build_dat(ntest = 20)
  res <- decal(dat$pert, dat$count, dat$clone)
  ## Generic checks
  check_decal(dat, res)
  ## Significance tests
  expect_equal(mean(res$pvalue < 0.05), 0.05, tolerance = 2 * 1 / 20)
  expect_equal(mean(res$p_adjusted < 0.05), 0, tolerance = 2 * 1 / 20)
})

test_that("decal recognize named indexes", {
  dat <- build_dat(named = TRUE)
  check_decal(dat, decal(dat$pert, dat$count, dat$clone))
})

test_that("decal accepts integer index with named matrix", {
  dat <- build_dat()

  cname <- sprintf("cell-%03d", seq_len(ncol(dat$count)))
  rname <- sprintf("gene-%03d", seq_len(nrow(dat$count)))
  gname <- sprintf("clone-%03d", seq_along(dat$cclone))
  dimnames(dat$count) <- list(rname, cname)

  check_decal(dat, decal(dat$pert, dat$count, dat$clone))
})

test_that("decal validates perturbation format", {
  count <- matrix(rpois(100 * 100, 10), ncol = 100)
  clone <- split(seq_len(100), sample(1:5, 100, replace = TRUE))
  perturbations <- data.frame(
    x = sample(which(rowMeans(count) > 1), 5, replace = TRUE),
    y = sample(1:5, 5, replace = TRUE)
  )
  expect_error(decal(integer(), count, clone), "must be a data.frame")
  expect_error(decal(matrix(), count, clone), "must be a data.frame")
  expect_error(
    decal(perturbations, count, clone),
    "must contain all required columns"
  )
})

test_that("decal validates count as matrix", {
  count <- matrix(rpois(100 * 100, 10), ncol = 100)
  clone <- split(seq_len(100), sample(1:5, 100, replace = TRUE))
  perturbations <- data.frame(
    gene = sample(which(rowMeans(count) > 1), 5, replace = TRUE),
    clone = sample(1:5, 5, replace = TRUE)
  )
  expect_error(
    decal(perturbations, integer(), clone),
    "must be a matrix"
  )
  expect_error(
    decal(perturbations, data.frame(), clone),
    "must be a matrix"
  )
})

test_that("decal validates clone as list", {
  count <- matrix(rpois(100 * 100, 10), ncol = 100)
  clone <- split(seq_len(100), sample(1:5, 100, replace = TRUE))
  perturbations <- data.frame(
    gene = sample(which(rowMeans(count) > 1), 5, replace = TRUE),
    clone = sample(1:5, 5, replace = TRUE)
  )
  expect_error(
    decal(perturbations, count, integer()),
    "must be a list"
  )
  expect_error(
    decal(perturbations, count, matrix()),
    "must be a list"
  )
})

test_that("decal output is consistent with perturbations", {
  count <- matrix(rpois(100 * 100, 10), ncol = 100)
  clone <- split(seq_len(100), sample(1:5, 100, replace = TRUE))
  perturbations <- data.frame(
    x = sample(which(rowMeans(count) > 1), 100, replace = TRUE),
    y = sample(1:5, 10, replace = TRUE)
  )
  expected_cols <- c(
    "x", "y", "n0", "n1", "x0", "x1", "mu", "theta",
    "xb", "z", "lfc", "pvalue", "p_adjusted"
  )

  result <- decal(perturbations, count, clone, gene_col = "x", clone_col = "y")
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
  expected_cols <- c(
    "n0", "n1", "x0", "x1", "mu", "theta",
    "xb", "z", "lfc", "pvalue", "p_adjusted"
  )
  ## forcing tests to fail
  expect_warning(
    fit <- decal(d, Y, X, min_n = 1E3),
    "No gene & clone perturbation matched the requirements"
  )
  expect_type(fit, typeof(d))
  expect_named(fit, c(names(d), expected_cols), ignore.order = TRUE)
  expect_equal(sum(is.na(fit$pvalue)), 10)
})
