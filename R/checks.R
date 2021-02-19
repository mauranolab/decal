#' @importFrom methods is
#' @noRd
is_matrix <- function(x) is.matrix(x) || is(x, "sparseMatrix")

#' @noRd
validate_numeric <- function(x) {
  if (!is.numeric(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a numeric", call. = FALSE)
  }
}

validate_numeric_vector <- function(x) {
  if (!is.numeric(x) || is.matrix(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a numeric vector", call. = FALSE)
  }
}

#' @noRd
validate_numeric_scalar <- function(x) {
  if (!is.numeric(x) || length(x) != 1) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a numeric scalar", call. = FALSE)
  }
}

#' @noRd
validate_integer_scalar <- function(x) {
  if (length(x) != 1 || !(is.numeric(x) && x == round(x))) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a integer scalar", call. = FALSE)
  }
}

#' @noRd
validate_positive_integer_scalar <- function(x) {
  if (length(x) != 1 || !(is.numeric(x) && x == round(x) && x > 0)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a positive integer scalar", call. = FALSE)
  }
}
