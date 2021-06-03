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
#' @param ... currently ignored parameters
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
#' @importFrom Matrix colSums rowMeans
#' @export
decal <- function(
  perturbations, count, clone, ..., theta_sample = 2000,
  min_mu = 0.05, min_n = 3, min_x = 1,
  gene_col="gene", clone_col="clone"
) {
  ## Validate input
  ## Build clone identity matrix
  cellnames <- colnames(count)
  if (is.null(cellnames)) { cellnames <- seq_len(ncol(count)) }
  clonemtx <- build_clone_matrix(clone, cellnames)
  ## Compute metrics
  n1 <- colSums(clonemtx)
  n0 <- ncol(count) - n1
  log_depth <- log(colSums(count))
  mu <- rowMeans(count)
  xp <- count %*% clone
  x1 <- coldiv(xp, n1)
  x0 <- coldiv(rowSums(count) - xp, ncol(count) - n1)
  theta <- estimate_theta(
    count, mu, log_depth, n = theta_sample, genes = which(mu >= min_mu)
  )
  ## Extract indexes
  rowidx <- get_index(perturbations[[gene_col]], rownames(count))
  colidx <- get_index(perturbations[[clone_col]], colnames(clone))
  mtxidx <- cbind(rowidx, colidx)
  if(any(is.na(rowidx)) || any(is.na(colidx))) {
    warning(
      "Skipping ", sum(is.na(rowidx) | is.na(colidx)), " perturbations ",
      "not matching gene or clone id", call.=FALSE)
  }
  ## Update results
  perturbations <- cbind(perturbations, data.frame(
    n0 = n0[colidx],
    n1 = n1[colidx],
    x1 = x1[mtxidx],
    x0 = x0[mtxidx],
    mu = mu[rowidx],
    theta = mu[rowidx],
    xb = NA_real_, z = NA_real_, lfc = NA_real_,
    pvalue = NA_real_, p_adjusted = NA_real_
  ))
  ## Fit analysis
  which_test <- which(
    (perturbations$x1 >= min_x | perturbations$x0 >= min_x) &
      perturbations$n1 >= min_n &
      perturbations$mu >= min_mu
  )
  if (length(which_test) == 0) {
    warning("No gene & clone perturbation matched the required criteria.")
  } else {
    fit <- fit_nb(count, clone, theta, log_depth,
                  rowidx[which_test], colidx[which_test])
    perturbations[which_test, names(fit)] <- fit
  }
  return(perturbations)
}

coldiv <- function(mtx, vec) {
  ## Validate input size
  if (length(vec) != ncol(mtx)) {
    stop("Incompatible dimensions. `vec` must have the same length as `mtx` columns",
         call. = FALSE)
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

name_or_index <- function(x) {
  name <- names(x)
  if (is.null(name)) { return(seq_along(x)) }
  return(name)
}

is_index <- function(x) {
  return(is.character(x) || (is.integer(x) && all(x > 0)))
}

get_index <- function(x, ref=NULL) {
  if (!is.null(ref) && typeof(x) == typeof(ref)) return(match(x, ref))
  if (is.integer(x) && all(x > 0L)) return(x)
  stop("`x` must be an integer index or match reference type")
}

validate_index <- function(x, unique=TRUE) {
  if (!is_index(x)) {
    stop("An index must be a character or integer vector", call. = FALSE)
  }
  if (unique && any(table(x) > 1L)) {
    stop("An index must not contain repeated entries", call. = FALSE)
  }
  return()
}

#' Build a clone identity matrix from a list of cells
#'
#' @importFrom Matrix sparseMatrix
#' @noRd
build_clone_matrix <- function(clone, cells) {
  ## Validate input
  if (!is.list(clone)) { stop("`clone` must be a list of cells", call. = FALSE) }
  validate_index(cells)
  ## Vectorize clone list
  clones <- name_or_index(clone)
  vec_cells <- unlist(clone)
  vec_clone <- rep(clones, sapply(clone, length))
  ## Build matrix
  mtx <- sparseMatrix(
    i = match(vec_cells, cells), j = match(vec_cells, clones), x=1L,
    dims=c(length(cells), length(clones))
  )
  if (is.character(cells)) rownames(mtx) <- cells
  if (is.character(clones)) colnames(mtx) <- clones
  return(mtx)
}

#' @importFrom fastglm fastglm
#' @importFrom MASS theta.ml
#' @importFrom stats poisson
#'
#' @noRd
estimate_theta_raw <- function(count, genes, log_dp) {
  fit_poisson <- function(i) {
    y <- Y[i,]
    x <- matrix(1L, nrow=length(y))
    f <- fastglm(x, y, family=poisson(), method=2, offset=log_dp)
    as.numeric(theta.ml(y, f$fitted.values))
  }
  Y <- as.matrix(count[genes,])

  theta <- suppressWarnings(sapply(seq_along(genes), FUN = fit_poisson))
  theta <- pmax(theta, 1E-7)
}

#' @importFrom stats approx bw.SJ density ksmooth predict
#'
#' @noRd
estimate_theta <- function(count, mu, log_dp, n, genes) {
  ## Response vector
  theta <- numeric(nrow(count))
  ## Subsample
  if (is.null(genes)) genes <- which(mu > 0)
  log_mu <- log(mu)[genes]
  ## Step1. Estimate theta in subsample
  ## TODO: subsample cols
  ## TODO: run estimate_theta_raw in blocks to avoid memory overflow
  if (n < length(genes)) {
    dens <- density(x=log_mu, bw="nrd", adjust=1)
    prob <- 1 / (approx(dens$x, dens$y, xout=log_mu)$y + .Machine$double.eps)
    gene1 <- sample(genes, size = n, prob = prob)
    log_mu1 <- log(mu)[gene1]
  } else {
    gene1 <- genes
    log_mu1 <- log_mu
  }
  theta1 <- estimate_theta_raw(count, genes, log_dp)
  ## Step2. Regularize theta estimate
  ## TODO: remove outliers?
  odds1 <- log10(1+10^log_mu1 /theta1)
  xs <- pmin(pmax(log_mu, min(log_mu1)), max(log_mu1))
  odds <- numeric(length(log_mu))
  odds[order(xs)] <- ksmooth(
    x = log_mu1, y = odds1, x.points = xs,
    bandwidth = bw.SJ(log_mu1) * 3, kernel = "normal"
  )$y
  theta[genes] <- 10^xs / (10^odds - 1)
  return(theta)
}

#' @importFrom fastglm fastglm
#' @importFrom MASS negative.binomial
#' @importFrom stats p.adjust
#'
#' @noRd
fit_nb <- function(count, clone, theta, log_dp, rows, cols) {
  if (length(theta) == 1L) theta <- rep(theta, nrow(Y))
  if (length(rows) != length(cols)) {
    stop("rows and cols index must have the same length", call. = FALSE)
  }
  ## Log-scale correction
  ln2 <- log(2)
  ## Subset to a dense matrix for speed
  ## TODO: blockade this step
  ix <- sort(unique(rows))
  jx <- sort(unique(cols))
  X <- as.matrix(clone[, jx, drop=FALSE])
  Y <- as.matrix(count[ix, , drop=FALSE])
  mean_depth <- mean(exp(log_dp))
  ## fit regression
  result <- mapply(function(i, j) {
    f <- fastglm(cbind(1, X[,j]), Y[i,],
                 family=negative.binomial(theta = theta[i]),
                 method=2, offset = log_dp)
    coef <- summary(f)$coef[3,]
    names(coef) <- NULL
    c(
      xb = predict(f, cbind(1, 1), type="response") * mean_depth,
      z = coef[3], lfc = coef[1]/ln2, pvalue = coef[4]
    )
  }, match(rows, ix), match(cols, jx))
  result <- as.data.frame(t(result))
  ## Compute Benjamin-Hochenberg p-value adjustment
  result$p_adjusted <- p.adjust(result$pvalue, "BH")
  return(result)
}
