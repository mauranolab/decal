#' Simulate a random clonal assignment
#'
#' It generates a random clone assignment list of specific sizes.
#'
#' @return a list of cells per clone.
#' @name simulate_clone
NULL

#' @rdname simulate_clone
#'
#' @param size number of cells to be assigned to each clone
#' @param ncells total number of cells available  to assign
#'
#' @importFrom stats setNames
#' @export
sim_clone <- function(size, ncells) {
  validate_positive_integer(size)
  validate_positive_integer_scalar(ncells)
  if (sum(size) > ncells) {
    stop("The total clones size exceed the total number of cells",
         call. = FALSE)
  }
  ci <- rep(seq_along(size), size)
  ri <- sample(seq_len(ncells), length(ci))
  return(setNames(split(ri, ci), NULL))
}

#' @rdname simulate_clone
#'
#' @param n total number of clones be generated
#' @param ncells total number of cells available to assign
#' @param min minimum number of cells perturbed by clone
#' @param max maximum number of cells perturbed by clone
#'
#' @export
sim_clone_range <- function(n, ncells, min = 2L, max = 20L) {
  validate_positive_integer_scalar(n)
  validate_positive_integer_scalar(ncells)
  validate_positive_integer_scalar(min)
  validate_positive_integer_scalar(max)
  if (min > max) {
    stop("Minimum number of cells per clone must be less than the maximum", call. = FALSE)
  }
  if (n * min > ncells) {
    stop("Unable to generate clones with current `ncells` limit", call. = FALSE)
  }

  size <- integer(n)
  for (i in seq(0, n - 1)) {
    limit <- pmin(max, (ncells - sum(size)) / (n-i))
    size[i + 1] <- sample(min:limit, 1)
  }
  sim_clone(size, ncells)
}
