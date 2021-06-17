#' Simulate a random UMI count matrix
#'
#' It generates a random UMI count matrix based on **DECAL** model assumptions
#' given the cells total count (`depth`) and gene prevalence ratio (ratio of
#' genes UMIs over the total; `ratio`).
#' Alternatively, you can further specify the generative model, by giving a
#' log2 fold-change (`lfc`) effect and each gene dispersion (`theta`).
#'
#' @return random UMI count matrix
#' @name simulate_count
NULL

#' @rdname simulate_count
#'
#' @param depth simulated cells total UMI count
#' @param ratio ratio of genes UMI to be simulated
#' @param lfc log2 fold-change effect to be applied to the simulated cells
#' @param theta genes' dispersion parameter
#'
#' @importFrom stats rnbinom
#' @export
sim_count <- function(depth, ratio, lfc = 0L, theta = 100L) {
  validate_positive_integer(depth)
  validate_rate(ratio)
  validate_numeric(lfc)
  validate_numeric(theta)
  ## Validate format
  if (length(theta) != 1L && length(theta) != length(ratio)) {
    stop("`theta` must be a scalar or vector of same length as `ratio`",
         call. = FALSE)
  }
  ## compute expected count
  if (is_matrix(lfc)) {
    if (nrow(lfc) != length(ratio) || ncol(lfc) != length(depth)) {
      stop("Incompatible dimensions. A `lfc` matrix must have the same",
           "number of rows as `ratio` length and columns as `depth` length",
           call. = FALSE)
    }
    mu <- (ratio %o% depth) * (2**lfc)
  } else if (length(lfc) == length(depth) || length(lfc) == 1L) {
    mu <- t((depth %o% ratio) * (2**lfc))
  } else {
    stop("`lfc` must be a scalar, vector of same length as `depth` or a matrix",
         call. = FALSE)
  }
  return(matrix(
    rnbinom(length(mu), size = theta, mu = mu),
    nrow = length(ratio), ncol = length(depth),
    dimnames = list(names(ratio), names(depth))
  ))
}

#' @rdname simulate_count
#'
#' @param reference a count matrix to base the simulation
#' @param lfc log2 fold-change effect to be applied to the simulated cells
#' @param theta genes' dispersion parameter
#' @param ngenes number of genes to be simulated. When `NULL` it replicates
#' the gene ratio found in `reference`, otherwise it simulates `ngenes`
#' with ratio ranging in a logarithm scale from lowest to highest observed
#' ratio observed.

#' @importFrom Matrix colSums rowSums
#' @export
sim_count_from_data <- function(reference, lfc = 0, theta = 100L, ngenes = NULL) {
  validate_matrix(reference)
  validate_numeric(lfc)
  validate_numeric(theta)

  depth <- colSums(reference)
  ratio <- rowSums(reference) / sum(depth)
  if (!is.null(ngenes)) {
    validate_positive_integer_scalar(ngenes)
    ratio <- seq_log(min(ratio), max(ratio), length_out = ngenes)
  }
  sim_count(depth, ratio, lfc = lfc, theta = theta)
}
