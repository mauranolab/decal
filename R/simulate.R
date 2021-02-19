#' Simulates a single-cell experiment based on our model
#'
#' It simulates an experiment based on our model assumptions
#'
#' @export
sim_expression <- function(cell_dp, gene_ps, log2_fc = 0L, theta = 100L) {
  validate_numeric_vector(cell_dp)
  validate_numeric_vector(gene_ps)
  validate_numeric(log2_fc)
  validate_numeric_vector(theta)

  if (is.matrix(log2_fc)) {
    if (nrow(log2_fc) != length(gene_ps) || ncol(log2_fc) != length(cell_dp)) {
      stop(
        "Incompatible dimensions. When `log2_fc` is a matrix, it must have ",
        "the same number of rows as `len(gene_ps)` and columns as ",
        "`len(cell_dp)`", call. = FALSE)
    }
    sim_mu <- (gene_ps %o% cell_dp) * (2 ** log2_fc)
  } else if (length(log2_fc) == length(cell_dp) || length(log2_fc) == 1L) {
    sim_mu <- t((cell_dp %o% gene_ps) * (2 ** log2_fc))
  } else {
    stop(
      "Incompatible dimensions. `log2_fc` must be a scalar or vector of ",
      "the same length as `cell_dp`", call. = FALSE)
  }
  sim_nb <- theta / (theta + sim_mu)
  return(matrix(
    rnbinom(length(sim_mu), size = theta, prob = sim_nb),
    nrow = length(gene_ps), ncol = length(cell_dp),
    dimnames = list(names(gene_ps), names(cell_dp))
  ))
}

#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @export
sim_from_data <- function(Y, log2_fc = 0L, theta = 100L) {
  if (!is_matrix(Y)) { stop("`Y` must be a matrix or sparseMatrix") }
  dp <- colSums2(Y)
  ps <- rowSums2(Y) / sum(dp)
  sim_expression(dp, ps, log2_fc, theta = theta)
}

#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @export
sim_from_data_as_exp_range <-
  function(Y, log2_fc = 0L, theta = 100L, n_genes = 100L) {
    if (!is_matrix(Y)) { stop("`Y` must be a matrix or sparseMatrix") }
    dp <- colSums2(Y)
    ps <- rowSums2(Y) / sum(dp)
    sim_ps <- seq_log(min(ps), max(ps), length_out = n_genes)
    sim_effect(dp, sim_ps, log2_fc, theta = theta)
  }

#' @export
sim_lfc_matrix <- function(log2_fc, clone_size, ncell, rep_each = 1L) {
  ## validation
  if (!is.numeric(log2_fc)) { stop("`log2_fc` must be a numeric value") }
  if (!is.integer(clone_size) && any(clone_size > max(ncell))) {
    stop("`clone_size` must be a integer value below or equal to ncell")
  }
  if (!is.integer(ncell) && length(ncell) != 1L && ncell <= 0) {
    stop("`ncell` must be a integer scalar bigger than 0")
  }
  if (!is.integer(rep_each) && all(rep_each <= 0)) {
    stop("`rep_each` must be positive integer bigger than 0")
  }
  if (length(log2_fc) != length(clone_size) && length(clone_size) != 1L) {
    stop("Incompatible sizes, `log2_fc` and `clone_size` must have the same length")
  }
  if (length(log2_fc) != length(rep_each) && length(rep_each) != 1L) {
    stop("Incompatible sizes, `log2_fc` and `rep_each` must have the same length")
  }
  ## fix values for compatibility
  if (length(clone_size) == 1L) {
    clone_size <- rep(clone_size, length(log2_fc))
  }
  if (length(rep_each) == 1L) {
    rep_each <- rep(rep_each, length(log2_fc))
  }
  do.call(rbind, lapply(seq_along(log2_fc), function(i) {
    vec <- sample(rep(c(0, log2_fc[i]), c(n - clone_size[i], clone_size[i])))
    t(matrix(vec , length(vec) , rep_each[i]))
  }))
}
