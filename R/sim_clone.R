#' Simulate a random clone matrix
#'
#' Given a vector of clone sizes (`perturbed`) and the total number of cells
#' (`n_cells`) and generates a random adjacency matrix with `n_cells` rows and
#' `length(perturbed)` columns.
#'
#' @param n_cells total number of cells
#' @return adjacency matrix indicating which cells (rows) belong to each clone
#' (column).
#'
#' @name simulate_clone
NULL

#' @rdname simulate_clone
#'
#' @param size number of cells to be assigned to each clone
#'
#' @export
sim_clone <- function(size, n_cells) {
  validate_positive_integer_vector(size)
  validate_positive_integer_scalar(n_cells)
  if (sum(size) > n_cells) {
    stop("The sum of all clones size must be fewer or equal to `n_cells`")
  }
  ci <- rep(seq_along(size), size)
  ri <- sample(seq_len(n_cells), length(ci))
  clones <- matrix(0, nrow = n_cells, ncol = length(size))
  clones[cbind(ri, ci)] <- 1
  return(clones)
}

#' @rdname simulate_clone
#'
#' @param n_clones total number of clones be generated
#' @param min minimum number of cells perturbed by clone
#' @param max maximum number of cells perturbed by clone
#'
#' @export
sim_clone_range <- function(n_clones, n_cells, min = 2L, max = 20L) {
  validate_positive_integer_scalar(n_clones)
  validate_positive_integer_scalar(n_cells)
  validate_positive_integer_scalar(min)
  validate_positive_integer_scalar(max)

  if (n_clones * min > n_cells) {
    stop("Unable to sample clone size with current `n_cells` limit",
      call. = FALSE
    )
  }

  size <- integer(n_clones)
  for (i in seq(0, n_clones - 1)) {
    limit <- pmin(max, (n_cells - sum(size)) / (n_clones - i))
    size[i + 1] <- sample(min:limit, 1)
  }
  sim_clone(size, n_cells)
}
