#' @noRd
is_integer <- function(x) is.integer(x) || (is.numeric(x) && all.equal(x, round(x)))
#' @noRd
is_index <- function(x) is.integer(x) || (is.numeric(x) && all(x > 0) && all.equal(x, round(x)))

#' @noRd
as_index <- function(x, ref = NULL) {
  if (is_index(x) && (is.null(ref) || !is_integer(ref))) {
    x
  } else if (typeof(x) == typeof(ref)) {
    match(x, ref)
  } else {
    stop("`x` must either be an integer or match `ref` type")
    NULL
  }
}

#' @noRd
col_div <- function(mtx, vec) {
  if (length(vec) != ncol(mtx)) {
    stop("Incompatible dimensions at `col_div`. `vec` must have the same length as `mtx` columns.")
  }
  if (inherits(x = mtx, what = "dgCMatrix")) {
    mtx@x <- mtx@x / vec[rep(seq_len(mtx@Dim[2]), diff(mtx@p))]
  } else if (inherits(x = mtx, what = "dgTMatrix")) {
    mtx@x <- mtx@x / vec[mtx@j + 1]
  } else {
    mtx <- t(t(mtx) / vec)
  }
  return(mtx)
}

#' Based on vst::get_model_pars function
#' (see: https://github.com/ChristophH/sctransform/blob/master/R/vst.R)
#'
#' @importFrom fastglm fastglm
#' @importFrom MASS theta.ml
#' @importFrom stats approx bw.SJ density ksmooth poisson predict
#' @noRd
estimate_theta <- function(count, mu, depth, n = 2e3, genes = NULL) {
  theta <- numeric(nrow(count))
  if (is.null(genes)) {
    genes <- which(mu > 0)
  }
  log_mean <- log(mu)[genes]
  ## Sample n genes to estimate
  genes_step1 <- genes
  log_mean_step1 <- log_mean
  if (n < length(genes_step1)) {
    log_mean_dens <- density(x = log_mean_step1, bw = "nrd", adjust = 1)
    sampling_prob <- 1 / (
      approx(log_mean_dens$x, log_mean_dens$y, xout = log_mean_step1)$y +
        .Machine$double.eps
    )
    genes_step1 <- sample(genes_step1, size = n, prob = sampling_prob)
    log_mean_step1 <- log(mu)[genes_step1]
  }
  ## Estimate theta
  Y <- as.matrix(count[genes_step1, ])
  theta_step1 <- suppressWarnings(sapply(seq_along(genes_step1), FUN = function(i) {
    y <- Y[i, ]
    fit <- fastglm(cbind(1, depth), y, family = poisson(), method = 2)
    as.numeric(theta.ml(y, fit$fitted.values))
  }))
  theta_step1 <- pmax(theta_step1, 1e-7)
  ## TODO: add find and remove outliers
  ## Regularize and predict
  odfac_step1 <- log10(1 + 10**log_mean_step1 / theta_step1)
  bw <- bw.SJ(log_mean_step1) * 3
  x_points <- pmin(pmax(log_mean, min(log_mean_step1)), max(log_mean_step1))
  o_points <- order(x_points)
  odfac <- numeric(length(log_mean))
  odfac[o_points] <- ksmooth(
    x = log_mean_step1, y = odfac_step1,
    x.points = x_points, bandwidth = bw, kernel = "normal"
  )$y
  theta[genes] <- 10**x_points / (10**odfac - 1)
  return(theta)
}

#' @importFrom fastglm fastglm
#' @importFrom MASS negative.binomial
#' @importFrom qvalue qvalue
#'
#' @noRd
fit_nb <- function(Y, X, theta, depth, is, js) {
  LN2 <- log(2) ## CTE
  if (length(theta) == 1) theta <- rep(theta, nrow(Y))
  ## Reindex to save memory
  ix <- sort(unique(is))
  jx <- sort(unique(js))
  X <- as.matrix(X[, jx, drop = FALSE])
  Y <- as.matrix(Y[ix, , drop = FALSE])
  mean_depth <- log(mean(exp(depth)))

  result <- as.data.frame(t(
    mapply(match(is, ix), match(js, jx), FUN = function(i, j) {
      fit <- fastglm(cbind(1, depth, X[, j]), Y[i, ],
        family = negative.binomial(theta = theta[i]), method = 2
      )
      coef <- summary(fit)$coef[3, ]
      names(coef) <- NULL
      c(
        xb = predict(fit, cbind(1, mean_depth, 1), type = "response"),
        beta = coef[1], stderr = coef[2], z = coef[3], lfc = coef[1] / LN2,
        pvalue = coef[4]
      )
    })
  ))
  ## TODO: consider move to BF as default and use `p.adjust`
  result$qvalue <- if (nrow(result) == 1) {
    result$pvalue
  } else if (nrow(result) < 1e3) {
    qvalue(result$pvalue, pi0 = 1)$qvalues
  } else {
    qvalue(result$pvalue)$qvalues
  }
  result
}

#' @noRd
validate_column <- function(x, name, ref) {
  if (!is_integer(x) && !is.character(x)) {
    stop("`", name, "` must be either a column name or index")
  } else if (is_integer(x) && (x > length(ref) || x <= 0)) {
    stop("`", name, "` is out of bounds")
  } else if (!is_integer(x) && !x %in% ref) {
    stop(
      "`", x, "` is not a column name. Set `",
      name, "` to use a different column"
    )
  }
}

#' Perturbed clone differential expression
#'
#' It takes scRNA-seq UMI count matrix, clone assignment matrix and compares
#' the expression of a given gene within and outside a given clone.
#'
#' The UMI counts are modeled by a Negative Binomial regression
#' (\eqn{Y_ij ~ NB(\mu_ij, \theta_i)}), where \eqn{\theta_i} is the gene \eqn{i}
#' dispersion and \eqn{\mu_ij} is gene \eqn{i} and cell \eqn{j} expected UMI
#' count estimated as follow:
#' \deqn{log(\mu_ij) = \beta_0 + \beta_1 \times log(dp_j) + \beta_2 x_ij}
#' where \eqn{dp_j} is cell \eqn{j} total UMI depth and \eqn{x_ij} is a
#' \emph{i.i.d.} measure if cell \eqn{j} belonged to a clone perturbing gene
#' \eqn{i} expression.
#'
#' @param pertubed A data.frame like structure describing clone and gene pairs
#'   to model expression and measure perturbation.
#' @param count A `j` x `z` logical matrix or sparseMatrix indicating if the
#'   \emph{jth} cell belongs or not to \emph{zth} clone.
#' @param clone A `i` x `j` count matrix or sparseMatrix indicating if the UMI
#'   expression of \emph{ith} gene for the  \emph{jth} cell.
#' @param min_x Minimal average expression within or outside the pertubed clone
#'   (`x1` and `x0`, respectively) required to test effect.
#' @param min_n Minimal number of cells required to test a clone.
#' @param theta_min_mu Minimal `mu` (overall average expression) required to
#'   estimated gene `theta` and to be tested.
#' @param theta_n Number of genes sampled to estimate `theta`.
#' @param gene_col Name or index of `gene` column in `pertubed`.
#' @param clone_col Name or index of `clone` column in `pertubed`.
#' @param ... Currently ignored
#'
#' @return A object of the same type as `pertubed` added the following columns:
#'
#' @importFrom DelayedMatrixStats colSums2 rowMeans2 rowSums2
#' @export
clone_nb <- function(pertubed, count, clone, ...,
                     min_x = 1, min_n = 2, theta_min_mu = 0.05, theta_n = 2000,
                     gene_col = "gene", clone_col = "clone") {
  ## Validate perturbed
  validate_dataframe(pertubed)
  validate_column(gene_col, "gene_col", colnames(pertubed))
  validate_column(clone_col, "clone_col", colnames(pertubed))
  ## Validate clone and count
  validate_matrix(count)
  validate_matrix(clone)
  if (any(clone > 1)) {
    stop("`clone` must be a logical matrix or encoded as 0/1")
  }
  if (ncol(count) != nrow(clone)) {
    stop("`ncol(count) != nrow(clone)`: have incompatible dimensions")
  }
  ## Validate filters
  validate_numeric_scalar(min_x)
  validate_numeric_scalar(min_n)
  validate_numeric_scalar(theta_min_mu)
  validate_numeric_scalar(theta_n)
  ## Extract indexes
  is <- as_index(pertubed[[gene_col]], rownames(count))
  js <- as_index(pertubed[[clone_col]], colnames(clone))
  ij <- cbind(is, js)
  # Pre-compute some metrics
  N <- colSums2(clone)
  D <- colSums2(count)
  Mu <- rowMeans2(count)
  Nz <- rowSums2(count > 0)
  XP <- count %*% clone
  X1 <- col_div(XP, N)
  X0 <- col_div(rowSums2(count) - XP, ncol(count) - N)
  Th <- estimate_theta(
    count, Mu, log(D),
    n = floor(theta_n), genes = which(Mu >= theta_min_mu)
  )

  pertubed[["n"]] <- as.integer(N[js])
  pertubed[["nonzero"]] <- as.integer(Nz[is])
  pertubed[["x1"]] <- X1[ij]
  pertubed[["x0"]] <- X0[ij]
  pertubed[["mu"]] <- Mu[is]
  pertubed[["theta"]] <- Th[is]

  to_test <- (
    (pertubed[["x1"]] >= min_x | pertubed[["x0"]] >= min_x) &
      pertubed[["n"]] >= min_n &
      pertubed[["mu"]] >= theta_min_mu
  )
  if (sum(to_test) == 0) {
    warning("No gene+clone passed the filter criteria, thus none were evaluated")
    return(pertubed)
  }
  fit <- fit_nb(count, clone, Th, log(D), is[to_test], js[to_test])
  pertubed[, names(fit)] <- NA_real_
  pertubed[to_test, names(fit)] <- fit
  pertubed
}
