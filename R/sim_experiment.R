#' Simulate a clone perturbation experiment
#'
#' Given a reference single-cell UMI count matrix, it generates a random set of
#' clones, a random count matrix, a table of clone and gene interactions and
#' the applied log2 fold-change.
#'
#' First, it produces a clone matrix using `sim_clone_range` to produce clones
#' size between `min_n` and `max_n`. Next, for each clone it generates a gene
#' interaction for each `log2_fc`. Finally, it computes the intended log2
#' fold-change effect and generates a random count matrix.
#'
#' @param reference reference UMI count matrix
#' @param log2_fc   log2 fold-change to apply to each interaction of a clone
#' @param n_clones    number of clones to generate
#' @param n_genes     number of genes to generate, if `NULL` it will generate
#' the same number of genes as `reference`
#' @param min_n     minimum clone size
#' @param max_n     maximum clone size
#' @param theta     genes' negative binomial dispersion
#' @return list containing table of gene and clone `interactions`, `clone`
#' matrix and `count` matrix.
#'
#' @export
sim_experiment <- function(reference, log2_fc = 0L, n_clones = 10L, n_genes = NULL,
                           min_n = 2L, max_n = 20L, theta = 100L) {
  validate_matrix(reference)
  validate_numeric_vector(log2_fc)
  validate_positive_integer_scalar(n_clones)
  if (is.null(n_genes) || length(n_genes) == 0 ||
    (is.numeric(n_genes) && length(n_genes) == 1L && n_genes > nrow(reference))) {
    nrows <- nrow(reference)
  } else {
    validate_positive_integer_scalar(n_genes)
    nrows <- n_genes
  }
  ncols <- ncol(reference)
  n_int <- length(log2_fc)
  if (n_int > nrows) {
    stop(
      "number of genes to be sampled for each clone exceed the total",
      .call = FALSE
    )
  }
  # produce clone matrix
  clone <- sim_clone_range(n_clones, ncols, min = min_n, max = max_n)
  # create interactions pairs
  g_vec <- seq_len(nrows)
  interactions <- data.frame(
    clone = rep(seq_len(n_clones), each = n_int),
    gene = as.vector(replicate(n_clones, sample(g_vec, n_int))),
    log2fc = rep(log2_fc, times = n_clones)
  )
  # compute lfc
  lfc <- matrix(0, nrow = nrows, ncol = ncols)
  for (i in seq_len(nrow(interactions))) {
    ci <- interactions$clone[i]
    ri <- interactions$gene[i]
    fc <- interactions$log2fc[i]
    if (any(lfc[ri, clone[, ci] == 1] != 0)) {
      warning(
        "Overlapping clone + gene effect was found.",
        "The intended log2 fold-change may be affected"
      )
    }
    lfc[ri, ] <- lfc[ri, ] + clone[, ci] * fc
  }
  count <- sim_count_from_data(reference, lfc, theta = theta, n_genes = n_genes)
  return(list(interactions = interactions, clone = clone, count = count))
}
