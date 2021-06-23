#' DECAL: Differential Expression analysis of Clonal Alterations Local effects
#' based on Negative Binomial distribution
#'
#' This function performs the clonal alterations differential expression
#' analysis pairs of clonal sub-populations and perturbed genes.
#'
#' Given a table of clone and gene pairs, a UMI count matrix, and list of cells
#' per clone, this function models gene expression (`Y`) with a negative
#' binomial (a.k.a. Gamma-Poisson) distribution for each perturbation pair as
#' a function of `X` (clone indicator variable) offset by the cell total count
#' (`D`) as described by the model:
#'
#' \deqn{Y ~ NB(xb, theta)}
#' \deqn{log(xb) = \beta_0 + \beta_x * X + log(D)}
#' \deqn{theta ~ \mu}
#'
#' The gene dispersion parameter (`theta`) is estimated and regularized in two
#' steps as developed by Hafemeister & Satija (2019). First, for a subset of
#' genes it fits a _Poisson_ regression offseted by `log(D)` and estimate a
#' crude `theta` using a maximum likelihood estimator with the observed counts
#' and regression results. Next, it regularize and expands `theta` estimates
#' with a kernel smoothing function as a function of average count (`mu`).
#'
#' @param perturbations table with clone and gene perturbations pairs to model
#' differential expression effect.
#' @param count UMI count matrix with cells as columns and genes (or features)
#' as rows.
#' @param clone list of cells per clone.
#' @param theta_sample number of genes sampled to preliminary `theta` estimation.
#' @param min_mu minimal overall average expression (`mu`) required.
#' @param min_n minimal number of perturbed cells (`n1`) required.
#' @param min_x minimal average expression of perturbed (`x1`) and non-perturbed
#' cells (`x0`) required.
#' @param gene_col gene index column name in `perturbations`
#' @param clone_col clone index column name in `perturbations`
#' @param p_method p-value adjustment for multiple comparisons.
#' See `\link[stats]{p.adjust}`.
#' @return it extends `perturbations` table adding the following columns:
#' - `n0` and `n1`: number of non-perturbed and perturbed cells
#' - `x0` and `x1`: number of non-perturbed and perturbed cells average count
#' - `mu`: overall average expression
#' - `theta`: negative binomial dispersion parameter
#' - `xb`: perturbed cells' estimated average count
#' - `z`: perturbed cells' standardize z-score effect
#' - `lfc`: perturbed cells' log2 fold-change effect
#' - `pvalue`
#' - `p_adjusted`
#'
#' @importFrom Matrix colSums rowSums rowMeans
#' @export
decal <- function(perturbations, count, clone, theta_sample = 2000,
                  min_mu = 0.05, min_n = 3, min_x = 1,
                  gene_col = "gene", clone_col = "clone", p_method = "BH") {
  ## Validate input
  validate_dataframe(perturbations, c(gene_col, clone_col))
  validate_matrix(count)
  ## Build clone identity matrix
  cellnames <- or(colnames(count), seq_len(ncol(count)))
  clonemtx <- build_clone_matrix(clone, cellnames)
  ## Compute metrics
  n1 <- colSums(clonemtx)
  n0 <- ncol(count) - n1
  xp <- count %*% clonemtx
  x1 <- coldiv(xp, n1)
  x0 <- coldiv(rowSums(count) - xp, ncol(count) - n1)
  mu <- rowMeans(count)
  log_depth <- log(colSums(count))
  theta <- estimate_theta(count, mu, log_depth,
    n = theta_sample,
    genes = which(mu >= min_mu)
  )
  ## Extract indexes
  rowidx <- get_index(perturbations[[gene_col]], rownames(count))
  colidx <- get_index(perturbations[[clone_col]], names(clone))
  mtxidx <- cbind(rowidx, colidx)
  if (any(is.na(rowidx)) || any(is.na(colidx))) {
    warning(
      "Skipping ", sum(is.na(rowidx) | is.na(colidx)), " perturbations ",
      "not matching gene or clone id"
    )
  }
  ## Update results
  perturbations <- cbind(perturbations, data.frame(
    n0 = ncol(count) - n1[colidx], n1 = n1[colidx],
    x0 = x0[mtxidx], x1 = x1[mtxidx], mu = mu[rowidx],
    theta = theta[rowidx],
    xb = NA_real_, z = NA_real_, lfc = NA_real_,
    pvalue = NA_real_, p_adjusted = NA_real_
  ))
  attr(perturbations, "raw_theta") <- attr(theta, "raw")
  ## Fit analysis
  which_test <- which(
    (perturbations$x1 >= min_x | perturbations$x0 >= min_x) &
      perturbations$n1 >= min_n &
      perturbations$mu >= min_mu
  )
  if (length(which_test) == 0) {
    warning("No gene & clone perturbation matched the requirements")
  } else {
    fit <- fit_nb(
      count, clonemtx, theta, log_depth,
      rowidx[which_test], colidx[which_test],
      p_method
    )
    perturbations[which_test, names(fit)] <- fit
  }
  return(perturbations)
}

#' @importFrom Matrix sparseMatrix
#' @noRd
build_clone_matrix <- function(clone, cells) {
  ## Validate input
  if (!is.list(clone)) {
    stop("`clone` must be a list of cells", call. = FALSE)
  }
  if (any(!is_index(unlist(clone)))) {
    stop("`clone` cells list must be integer or character indexes", call. = FALSE)
  }
  validate_index(cells)
  ## Vectorize clone list
  clones <- or(names(clone), seq_along(clone))
  vec_cells <- unlist(clone)
  vec_clone <- rep(clones, sapply(clone, length))
  ## Build matrix
  mtx <- sparseMatrix(
    i = match(vec_cells, cells), j = match(vec_clone, clones), x = 1L,
    dims = c(length(cells), length(clones))
  )
  if (is.character(cells)) rownames(mtx) <- cells
  if (is.character(clones)) colnames(mtx) <- clones
  return(mtx)
}

#' @importFrom fastglm fastglm
#' @importFrom MASS theta.ml
#' @importFrom stats poisson
#' @noRd
estimate_theta_raw <- function(count, genes, log_dp, floor = 1E-7) {
  Y <- as.matrix(count[genes, ])
  theta <- sapply(seq_along(genes), FUN = function(i) {
    y <- Y[i, ]
    x <- matrix(1, nrow = length(y))
    f <- fastglm(x, y, family = poisson(), method = 2, offset = log_dp)
    suppressWarnings(as.numeric(theta.ml(y, f$fitted.values)))
  })
  return(pmax(theta, floor))
}

#' @importFrom stats bw.SJ ksmooth
#' @noRd
regularize_theta <- function(logmu1, theta1, logmu) {
  odds1 <- log10(1 + 10**logmu1 / theta1)
  xpnts <- pmin(pmax(logmu, min(logmu1)), max(logmu1))
  odds <- numeric(length(logmu))
  odds[order(xpnts)] <- ksmooth(
    x = logmu1, y = odds1, x.points = xpnts,
    bandwidth = bw.SJ(logmu1) * 3, kernel = "normal"
  )$y
  return(10**xpnts / (10**odds - 1))
}

#' @importFrom stats approx density
#' @noRd
estimate_theta <- function(count, mu, log_dp, n = 2000, genes) {
  ## Subsample
  genes <- or(genes, which(mu > 0))
  logmu <- log(mu)[genes]
  ## Step1. Estimate theta in subsample
  ## TODO: subsample cols
  ## TODO: run estimate_theta_raw in blocks to avoid memory overflow
  if (n < length(genes)) {
    dens <- density(x = logmu, bw = "nrd", adjust = 1)
    prob <- 1 / (approx(dens$x, dens$y, xout = logmu)$y + .Machine$double.eps)
    genes1 <- sample(genes, size = n, prob = prob)
    logmu1 <- log(mu)[genes1]
  } else {
    genes1 <- genes
    logmu1 <- logmu
  }
  theta1 <- estimate_theta_raw(count, genes, log_dp)
  ## Step2. Regularize theta estimate
  ## TODO: remove outliers?
  raw <- numeric(nrow(count))
  raw[genes1] <- theta1
  theta <- numeric(nrow(count))
  theta[genes] <- regularize_theta(logmu1, theta1, logmu)
  attr(theta, "raw") <- raw
  return(theta)
}

#' @importFrom fastglm fastglm
#' @importFrom MASS negative.binomial
#' @importFrom stats p.adjust predict
#' @noRd
fit_nb <- function(count, clone, theta, log_dp, rows, cols, p_method) {
  if (length(theta) == 1L) theta <- rep(theta, nrow(count))
  if (length(rows) != length(cols)) {
    stop("rows and cols index must have the same length", call. = FALSE)
  }
  ## Log-scale correction
  ln2 <- log(2)
  ## Subset to a dense matrix for speed TODO: blockade this step
  ix <- sort(unique(rows))
  jx <- sort(unique(cols))
  X <- as.matrix(clone[, jx, drop = FALSE])
  Y <- as.matrix(count[ix, , drop = FALSE])
  mean_depth <- mean(exp(log_dp))
  ## fit regression
  result <- mapply(function(i, j) {
    f <- fastglm(cbind(1, X[, j]), Y[i, ],
      family = negative.binomial(theta = theta[i]), method = 2,
      offset = log_dp
    )
    coef <- summary(f)$coef[2, ]
    names(coef) <- NULL
    c(
      xb = predict(f, cbind(1, 1), type = "response") * mean_depth,
      z = coef[3], lfc = coef[1] / ln2, pvalue = coef[4]
    )
  }, match(rows, ix), match(cols, jx))
  result <- as.data.frame(t(result))
  ## Compute Benjamin-Hochenberg p-value adjustment
  result$p_adjusted <- p.adjust(result$pvalue, p_method)
  return(result)
}
