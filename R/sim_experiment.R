#' Simulate a DECAL perturbation experiment
#'
#' It simulates a random DECAL experiment generating a random UMI count matrix
#' with cells randomly assigned to different clones and introducing to each
#' clone perturbations to random genes.
#'
#' @param lfc vector indicating the perturbations to be generated for each clone
#' @param nclones number of clones to be simulated
#' @param min_n minimum clone size
#' @param max_n maximum clone size
#' @param theta negative binomial dispersion parameter
#' @return a list containing the following fields:
#' - `perturbations`, a table indicating the perturbations `expected_lfc` for
#'    a `gene` in a `clone`.
#' - `clone`, a list of cells assigned to each clone.
#' - `count`, a random UMI count matrix.
#'
#' @name simulate_experiment
NULL

#' @rdname simulate_experiment
#'
#' @param depth cells total UMI count
#' @param ratio ratio of genes UMI over total UMI count
#' @export
sim_experiment <- function(depth, ratio, lfc=0, nclones=10L, min_n=2L, max_n=20L, theta = 100L) {
  validate_positive_integer(depth)
  validate_rate(ratio)
  validate_numeric(lfc)
  validate_positive_integer_scalar(nclones)
  validate_positive_integer_scalar(min_n)
  validate_positive_integer_scalar(max_n)
  validate_numeric(theta)

  ncells <- length(depth)
  ngenes <- length(ratio)
  ## Assign to clones
  clone <- sim_clone_range(nclones, ncells, min_n, max_n)
  ## Simulate random perturbations
  perturbations <- simulate_perturbations(nclones, ngenes, lfc)
  ## produce count matrix
  lfc_mat <- build_lfc_mat(perturbations, clone, ngenes, ncells)
  count <- sim_count(depth, ratio, lfc = lfc_mat, theta = theta)
  return(list(perturbations = perturbations, clone = clone, count = count))
}

#' @rdname simulate_experiment
#'
#' @param reference UMI count matrix used as base for the simulation
#' @param ngenes number of genes to be simulated. When `NULL` it replicates
#' the gene ratio found in `reference`, otherwise it simulates `ngenes`
#' with ratio ranging in a logarithm scale from lowest to highest observed
#' ratio observed.
#' @export
sim_experiment_from_data <- function(reference, lfc=0, nclones=10L, min_n=2L, max_n=20L, theta = 100L, ngenes = NULL) {
  validate_matrix(reference)
  validate_numeric(lfc)
  validate_positive_integer_scalar(nclones)
  validate_positive_integer_scalar(min_n)
  validate_positive_integer_scalar(max_n)
  validate_numeric(theta)

  ncells  <- ncol(reference)
  ngenes_ <- nrow(reference)
  if (!is.null(ngenes)) {
    validate_positive_integer_scalar(ngenes)
    ngenes_ <- ngenes
  }

  clone <- sim_clone_range(nclones, ncells, min_n, max_n)
  perturbations <- simulate_perturbations(nclones, ngenes_, lfc)
  ## Build count matrix
  lfc_mat <- build_lfc_mat(perturbations, clone, ngenes_, ncells)
  count <- sim_count_from_data(reference, lfc_mat, theta = theta, ngenes = ngenes)

  return(list(perturbations = perturbations, clone = clone, count = count))
}

#' @noRd
build_lfc_mat <- function(perturbations, clone, nrow, ncol) {
  clone_mat <- build_clone_matrix(clone, seq_len(ncol))
  lfc_mat <- matrix(0, nrow = nrow, ncol = ncol)
  for (i in seq_len(nrow(perturbations))) {
    clone_i <- perturbations$clone[i]
    gene_i  <- perturbations$gene[i]
    lfc_i   <- perturbations$expected_lfc[i]
    ## WARN PRE-ASSIGNED PERTURBATION WAS OVERWRITTEN DUE TWO CLONES
    ## PERTUBATING THE SAME GENE
    x <- clone_mat[, clone_i] == 1
    if (any(lfc_mat[gene_i, x] != 0)) {
      warning("Two or more clones pertubated the same gene.",
              "The intended log2 fold-change may be affected.")
    }
    lfc_mat[gene_i,] <- lfc_mat[gene_i,] + x * lfc_i
  }
  return(lfc_mat)
}

#' @noRd
simulate_perturbations <- function(nclones, ngenes, lfc) {
  genes <- seq_len(ngenes)
  npert <- length(lfc)
  if (npert > ngenes) {
    stop("number of perturbations exceed the total number of genes", call. = FALSE)
  }
  data.frame(
    clone = rep(seq_len(nclones), each = npert),
    gene  = as.vector(replicate(nclones, sample(genes, npert))),
    expected_lfc = rep(lfc, times = nclones)
  )
}
