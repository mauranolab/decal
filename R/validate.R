#' @noRd
is_index <- function(x) {
  return(is.character(x) || (all(x == as.integer(x)) && all(x > 0)))
}

#' @noRd
validate_index <- function(x, unique = TRUE) {
  if (!is_index(x)) {
    stop("An index must be a character or integer vector", call. = FALSE)
  }
  if (unique && any(table(x) > 1L)) {
    stop("An index must not contain repeated entries", call. = FALSE)
  }
}

#' @noRd
validate_dataframe <- function(x, cols = c()) {
  if (!is.data.frame(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a data.frame like structure (e.g. tibble)",
      call. = FALSE
    )
  }
  if (!all(cols %in% colnames(x))) {
    name <- deparse(substitute(x))
    stop("`", name, "` must contain all required columns: ",
      paste0("`", cols, "`", collapse = ", "), "\nfound:",
      paste0("`", colnames(x), "`", collapse = ", "),
      call. = FALSE
    )
  }
}

#' @importFrom methods is
#' @noRd
is_matrix <- function(x) {
  (!is.character(x) && is.matrix(x)) || is(x, "sparseMatrix")
}

#' @noRd
validate_matrix <- function(x) {
  if (!is_matrix(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a matrix or sparseMatrix", call. = FALSE)
  }
}

#' @noRd
is_integer <- function(x) {
  is.integer(x) || (is.numeric(x) && all(x == round(x)))
}

#' @noRd
validate_positive_integer <- function(x) {
  if (!(is_integer(x) && all(x > 0))) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a positive integer", call. = FALSE)
  }
}

#' @noRd
validate_positive_integer_scalar <- function(x) {
  if (!(is_integer(x) && all(x > 0) && length(x) == 1L)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a positive integer scalar", call. = FALSE)
  }
}

#' @noRd
validate_numeric <- function(x) {
  if (!is.numeric(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be numeric", call. = FALSE)
  }
}

#' @noRd
validate_numeric_scalar <- function(x) {
  if (!(is.numeric(x) && length(x) == 1L)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a numeric scalar", call. = FALSE)
  }
}

#' @noRd
validate_rate <- function(x) {
  do.call(validate_numeric, list(enquote(x)))
  if (!all(x >= 0 & x <= 1)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a ratio between 0 and 1", call. = FALSE)
  }
}
