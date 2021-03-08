#' Simulates a single-cell experiment UMI count matrix
#'
#' It simulates a random single-cell experiment UMI count matrix based on our
#' model assumptions given the total cell count (`cell_depth`) and gene rate
#' (genes total count divided by total count; `gene_rate`).
#' You can further specify the model, by giving a log2 fold-change effect to
#' apply by cell or by cell and gene (`log2_fc`) and each gene negative
#' binomial model dispersion (`theta`).
#'
#' @param log2_fc     log2 fold-change effect among the cells
#' @param theta       genes negative binomial dispersion parameter
#' @return random count matrix
#'
#' @name simulate_count
NULL

#' @rdname simulate_count
#'
#' @param cell_depth  total UMI count per cell
#' @param gene_rate   gene rate (gene total UMI divided by matrix total UMI)
#'
#' @importFrom stats rnbinom
#' @export
sim_count <- function(cell_depth, gene_rate, log2_fc = 0L, theta = 100L) {
  validate_numeric_vector(cell_depth)
  validate_numeric_vector(gene_rate)
  validate_numeric(log2_fc)
  validate_numeric_vector(theta)

  if (length(theta) != 1L && length(theta) != length(gene_rate)) {
    stop(
      "Incompatible dimensions. `theta` must be either a scalar or have the same
      length as `gene_rate`",
      call. = FALSE
    )
  }

  if (is_matrix(log2_fc)) {
    if (nrow(log2_fc) != length(gene_rate) || ncol(log2_fc) != length(cell_depth)) {
      stop(
        "Incompatible dimensions. When `log2_fc` is a matrix, it must have ",
        "the same number of rows as `len(gene_rate)` and columns as ",
        "`len(cell_depth)`",
        call. = FALSE
      )
    }
    sim_mu <- (gene_rate %o% cell_depth) * (2**log2_fc)
  } else if (length(log2_fc) == length(cell_depth) || length(log2_fc) == 1L) {
    sim_mu <- t((cell_depth %o% gene_rate) * (2**log2_fc))
  } else {
    stop(
      "Incompatible dimensions. `log2_fc` must be a scalar or vector of ",
      "the same length as `cell_depth`",
      call. = FALSE
    )
  }
  sim_nb <- theta / (theta + sim_mu)
  return(matrix(
    rnbinom(length(sim_mu), size = theta, prob = sim_nb),
    nrow = length(gene_rate), ncol = length(cell_depth),
    dimnames = list(names(gene_rate), names(cell_depth))
  ))
}

#' @rdname simulate_count
#'
#' @param reference   count matrix to be used as reference
#' @param n_genes     number of genes to simulate. If set to `NULL` it
#' produces the gene rate observed in `reference`, otherwise if a postive
#' integer scalar is given, it simulates `n_genes` with a gene rate spanning
#' from lowest to highest observed in a logarithm range observed.
#'
#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @export
sim_count_from_data <-
  function(reference, log2_fc = 0L, theta = 100L, n_genes = NULL) {
    validate_matrix(reference)
    dp <- colSums2(reference)
    ps <- rowSums2(reference) / sum(dp)
    if (!is.null(n_genes)) {
      validate_positive_integer_scalar(n_genes)
      ps <- seq_log(min(ps), max(ps), length_out = n_genes)
    }
    sim_count(dp, ps, log2_fc, theta = theta)
  }
