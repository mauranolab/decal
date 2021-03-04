#' Simulates a single-cell experiment based on our model
#'
#' It simulates an experiment based on our model assumptions
#'
#' @export
sim_count <- function(cell_dp, gene_ps, log2_fc = 0L, theta = 100L) {
  validate_numeric_vector(cell_dp)
  validate_numeric_vector(gene_ps)
  validate_numeric(log2_fc)
  validate_numeric_vector(theta)

  if (is_matrix(log2_fc)) {
    if (nrow(log2_fc) != length(gene_ps) || ncol(log2_fc) != length(cell_dp)) {
      stop(
        "Incompatible dimensions. When `log2_fc` is a matrix, it must have ",
        "the same number of rows as `len(gene_ps)` and columns as ",
        "`len(cell_dp)`",
        call. = FALSE
      )
    }
    sim_mu <- (gene_ps %o% cell_dp) * (2**log2_fc)
  } else if (length(log2_fc) == length(cell_dp) || length(log2_fc) == 1L) {
    sim_mu <- t((cell_dp %o% gene_ps) * (2**log2_fc))
  } else {
    stop(
      "Incompatible dimensions. `log2_fc` must be a scalar or vector of ",
      "the same length as `cell_dp`",
      call. = FALSE
    )
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
sim_count_from_data <-
  function(data, log2_fc = 0L, theta = 100L, n_genes = NULL) {
    validate_matrix(data)
    dp <- colSums2(data)
    ps <- rowSums2(data) / sum(dp)
    if (!is.null(n_genes)) {
      ps <- seq_log(min(ps), max(ps), length_out = n_genes)
    }
    sim_count(dp, ps, log2_fc, theta = theta)
  }

#' @export
sim_clone <- function(perturbed, cells, rep_each = NULL) {
  validate_positive_integer_vector(perturbed)
  validate_positive_integer_scalar(cells)
  if (!is.null(rep_each)) {
    validate_positive_integer_scalar(rep_each)
    perturbed <- rep(perturbed, rep_each)
  }
  if (any(perturbed > cells)) {
    stop("All `perturbed` values must smaller than `cells`")
  }
  sapply(seq_along(perturbed), function(i) {
    vec <- rep(0:1, c(cells - perturbed[i], perturbed[i]))
    sample(vec)
  })
}

#' #' @export
#' sim_pairs <- function(per_clone, genes, clones) {
#'   validate_positive_integer_scalar(per_clone)
#'   validate_positive_integer_scalar(genes)
#'   validate_positive_integer_scalar(clones)
#'   gene_pool <- seq_len(genes)
#'   data.frame(
#'     gene = replicate(length(clones), sample(gene_pool, per_clone))
#'     clone = rep(seq_len(clones), per_clone)
#'   )
#' }

#' @export
make_lfc <- function(clone, log2_fc = 1) {
  validate_adjacency_matrix(clone)
  validate_numeric_vector(log2_fc)
  ## validate dimensions
  do.call(rbind, lapply(seq_along(log2_fc), function(i) {
    t(clone * log2_fc[i])
  }))
}

# sim_lfc <- function(log2_fc, clone_size, ncell, rep_each = 1L) {
#   validate_numeric_vector(log2_fc)
#   validate_positive_integer_vector(clone_size)
#   validate_positive_integer_scalar(ncell)
#   validate_positive_integer_vector(rep_each)
#   if (any(clone_size > ncell)) stop("All `clone_size` values must smaller than `ncell`")
#   ## validate dimensions
#   len_lfc <- length(log2_fc)
#   if (length(clone_size) != 1L && length(clone_size) != len_lfc) {
#     stop(
#       "Incompatible dimensions. `clone_size` must be either a scalar ",
#       "or the same length as `log2_fc`"
#     )
#   }
#   if (length(rep_each) != 1L && length(rep_each) != len_lfc) {
#     stop(
#       "Incompatible dimensions. `rep_each` must be either a scalar ",
#       "or the same length as `log2_fc`"
#     )
#   }
#   ## Fix values for compatibility
#   if (length(clone_size) == 1L) clone_size <- rep(clone_size, len_lfc)
#   if (length(rep_each) == 1L) rep_each <- rep(rep_each, len_lfc)
#
#   do.call(rbind, lapply(seq_len(len_lfc), function(i) {
#     vec <- rep(c(0, log2_fc[i]), c(n - clone_size[i], clone_size[i]))
#     vec <- sample(vec) ## shuffle perturbed cells
#     t(matrix(vec, length(vec), rep_each[i]))
#   }))
# }
