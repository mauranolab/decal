#' @importFrom methods is
#' @noRd
is_matrix <- function(x) is.matrix(x) || is(x, "sparseMatrix")

#' @noRd
validate_matrix <- function(x) {
  if (!is_matrix(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a matrix or sparseMatrix", call. = FALSE)
  }
}

#' @importFrom methods is
#' @noRd
validate_adjacency_matrix <- function(x) {
  do.call(validate_matrix, list(enquote(x)))
  if (!is.logical(x) && any(x != 1L & x != 0L)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a logical matrix or encoded as 0/1",
      call. = FALSE
    )
  }
}

#' @noRd
validate_dataframe <- function(x) {
  if (!is.data.frame(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a data.frame like structure (e.g. tibble)",
      call. = FALSE
    )
  }
}

#' @noRd
validate_numeric <- function(x) {
  if (!is.numeric(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a numeric", call. = FALSE)
  }
}

validate_numeric_vector <- function(x) {
  do.call(validate_numeric, list(enquote(x)))
  if (is_matrix(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a numeric vector", call. = FALSE)
  }
}

#' @noRd
validate_numeric_scalar <- function(x) {
  do.call(validate_numeric, list(enquote(x)))
  if (length(x) != 1L) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a numeric scalar", call. = FALSE)
  }
}

#' @noRd
validate_integer_vector <- function(x) {
  if (!(is.numeric(x) && x == round(x)) || is_matrix(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a integer", call. = FALSE)
  }
}

#' @noRd
validate_positive_integer_vector <- function(x) {
  do.call(validate_integer_vector, list(enquote(x)))
  if (any(x <= 0)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a positive integer vector", call. = FALSE)
  }
}

#' @noRd
validate_integer_scalar <- function(x) {
  do.call(validate_integer_vector, list(enquote(x)))
  if (length(x) != 1) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a integer scalar", call. = FALSE)
  }
}

#' @noRd
validate_positive_integer_scalar <- function(x) {
  do.call(validate_integer_scalar, list(enquote(x)))
  if (x <= 0) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a positive integer scalar", call. = FALSE)
  }
}
