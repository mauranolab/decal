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
  if (is.null(genes)) genes <- which(mu > 0)
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
fit_nb <- function(Y, X, theta, log_depth, is, js) {
  if (length(theta) == 1) theta <- rep(theta, nrow(Y))
  ## Scale correction constant
  LN2 <- log(2)
  ## Subset and covert to dense matrix
  ix <- sort(unique(is))
  jx <- sort(unique(js))
  X <- as.matrix(X[, jx, drop = FALSE])
  Y <- as.matrix(Y[ix, , drop = FALSE])
  mean_depth <- log(mean(exp(log_depth)))

  result <- as.data.frame(t(
    mapply(function(i, j) {
      fit <- fastglm(
        cbind(1, log_depth, X[, j]), Y[i, ],
        family = negative.binomial(theta = theta[i]),
        method = 2
      )
      coef <- summary(fit)$coef[3, ]
      names(coef) <- NULL
      c(
        xb = predict(fit, cbind(1, mean_depth, 1), type = "response"),
        z = coef[3], lfc = coef[1] / LN2, pvalue = coef[4]
      )
    }, match(is, ix), match(js, jx))
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

#' Perturbed clone and gene differential expression
#'
#' Given a table of clone and gene interactions pairs, UMI count matrix and
#' clone adjacency matrix, it models each interaction gene expression (`Y`)
#' with a negative binomial regression in function of cell total count (`D`)
#' and a clone indicator variable (`X`) as indicated by the model below:
#'
#' \deqn{Y ~ NB(xb, theta)}
#' \deqn{log(xb) = \beta_0 + \beta_d * log(D) + \beta_x * X}
#' \deqn{theta ~ \mu}
#'
#' Each gene overdispersion (`theta`) is estimated in two steps. First, it fits
#' a Poisson regression as function of log total cell count and uses maximum
#' likelihood estimator to estimate the raw `theta` for a subset of genes.
#' Next it regularize and expands `theta` estimates as function of overall
#' average count.
#'
#' @param interactions table with clone and gene indexes pairs to measure
#' differential gene expression.
#' @param count count matrix with cells as columns and features as rows.
#' @param clone logical or adjacency matrix indicating the perturbed cells
#' (rows) belonging to each clone (columns).
#' @param min_x minimal average expression of perturbed and non-perturbed cells
#' required.
#' @param min_n minimal number of perturbed cells required.
#' @param theta_min_mu minimal overall average expression required to estimate
#' gene `theta` and run test.
#' @param theta_n number of genes sampled to preliminary `theta` estimation.
#' @param gene_col name or index of feature column in `interactions`
#' @param clone_col name or index of clone column in `interactions`
#' @param ... currently ignored
#' @return object of the type as `interactions` with the following columns added:
#'  - `n0` and `n1`: number of non-perturbed and perturbed cells
#'  - `x0` and `x1`: number of non-perturbed and perturbed cells average count
#'  - `mu`: overall average expression
#'  - `theta`: negative binomial overdispersion
#'  - `xb`: perturbed cells' estimated average count
#'  - `z`: perturbed cells' standardize z-score effect
#'  - `lfc`: perturbed cells' log2 fold-change effect
#'  - `pvalue`
#'  - `qvalue`
#'
#' @importFrom DelayedMatrixStats colSums2 rowMeans2 rowSums2
#' @export
clone_nb <- function(interactions, count, clone, ...,
                     min_x = 1, min_n = 2, theta_min_mu = 0.05, theta_n = 2000,
                     gene_col = "gene", clone_col = "clone") {
  ## Validate perturbed
  validate_dataframe(interactions)
  validate_column(gene_col, "gene_col", colnames(interactions))
  validate_column(clone_col, "clone_col", colnames(interactions))
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
  validate_positive_integer_scalar(min_n)
  validate_numeric_scalar(theta_min_mu)
  validate_positive_integer_scalar(theta_n)
  if (theta_n < 10L) {
    stop("`theta_n` should be larger than 10", .call = FALSE)
  }
  # Pre-compute some metrics
  N <- colSums2(clone)
  D <- log(colSums2(count))
  Mu <- rowMeans2(count)
  Nz <- rowSums2(count > 0)
  XP <- count %*% clone
  X1 <- col_div(XP, N)
  X0 <- col_div(rowSums2(count) - XP, ncol(count) - N)
  Th <- estimate_theta(
    count, Mu, D,
    n = theta_n,
    genes = which(Mu >= theta_min_mu)
  )
  ## Extract indexes
  is <- as_index(interactions[[gene_col]], rownames(count))
  js <- as_index(interactions[[clone_col]], colnames(clone))
  ij <- cbind(is, js)
  ## Add statistics
  interactions[["n0"]] <- as.integer(ncol(count) - N[js])
  interactions[["n1"]] <- as.integer(N[js])
  interactions[["x1"]] <- X1[ij]
  interactions[["x0"]] <- X0[ij]
  interactions[["mu"]] <- Mu[is]
  interactions[["theta"]] <- Th[is]
  interactions[, c("xb", "z", "lfc", "pvalue", "qvalue")] <- NA_real_

  which_to_test <- which(
    (interactions[["x1"]] >= min_x | interactions[["x0"]] >= min_x) &
      interactions[["n1"]] >= min_n &
      interactions[["mu"]] >= theta_min_mu
  )
  if (length(which_to_test) == 0) {
    warning("No gene+clone passed the filter criteria, thus none were evaluated")
    return(interactions)
  }
  fit <- fit_nb(count, clone, Th, D, is[which_to_test], js[which_to_test])
  interactions[which_to_test, names(fit)] <- fit
  interactions
}
