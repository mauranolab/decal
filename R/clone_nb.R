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
  if (inherits(x = mtx, what = 'dgCMatrix')) {
    mtx@x <- mtx@x / vec[rep(seq_len(mtx@Dim[2]), diff(mtx@p))]
  } else if (inherits(x = mtx, what = 'dgTMatrix')) {
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
#' @noRd
estimate_theta <- function(count, mu, depth, n=2e3, genes=NULL) {
  theta <- numeric(nrow(count))
  if(is.null(genes)) { genes <- which(mu > 0) }
  log_mean <- log(mu)[genes]
  ## Sample n genes to estimate
  genes_step1 <- genes
  log_mean_step1 <- log_mean
  if (n < length(genes_step1)) {
    log_mean_dens <- density(x=log_mean_step1, bw="nrd", adjust=1)
    sampling_prob <- 1 / (
      approx(log_mean_dens$x, log_mean_dens$y, xout=log_mean_step1)$y +
        .Machine$double.eps
      )
    genes_step1 <- sample(genes_step1, size=n, prob = sampling_prob)
    log_mean_step1 <- log(mu)[genes_step1]
  }
  ## Estimate theta
  Y <- as.matrix(count[genes_step1, ])
  theta_step1 <- suppressWarnings(sapply(seq_along(genes_step1), FUN = function(i) {
    y <- Y[i,]
    fit <- fastglm(cbind(1, depth), y, family=poisson(), method=2)
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
  odfac[o_points] <- ksmooth(x=log_mean_step1, y=odfac_step1,
    x.points=x_points, bandwidth=bw, kernel="normal")$y
  theta[genes] <- 10 ** x_points / (10 ** odfac - 1)
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
    mapply(match(is, ix), match(js, jx), FUN=function(i, j) {
      fit <- fastglm(cbind(1, depth, X[,j]), Y[i,],
        family=negative.binomial(theta = theta[i]), method=2)
      coef <- summary(fit)$coef[3,]
      names(coef) <- NULL
      c(
        xb = predict(fit, cbind(1, mean_depth, 1), type="response"),
        beta = coef[1], stderr = coef[2], z = coef[3], lfc = coef[1]/LN2,
        pvalue = coef[4]
      )
    })
  ))
  ## TODO: consider move to BF as default and use `p.adjust`
  result$qvalue <- if (nrow(result) == 1) {
    result$pvalue
  } else if (nrow(result) < 1e3) {
    qvalue(result$pvalue, pi0=1)$qvalues
  } else {
    qvalue(result$pvalue)$qvalues
  }
  result
}

