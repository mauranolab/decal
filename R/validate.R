#' @noRd
is_index <- function(x) {
  return(is.character(x) || (is.integer(x) && all(x > 0)))
}

#' @noRd
validate_index <- function(x, unique=TRUE) {
  if (!is_index(x)) {
    stop("An index must be a character or integer vector", call. = FALSE)
  }
  if (unique && any(table(x) > 1L)) {
    stop("An index must not contain repeated entries", call. = FALSE)
  }
  return()
}

#' @noRd
validate_dataframe <- function(x, cols=c()) {
  if (!is.data.frame(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a data.frame like structure (e.g. tibble)",
         call. = FALSE
    )
  }
  if (all(cols %in% colnames(x))) {
    name <- deparse(substitute(x))
    stop("`", name, "` must contain all required columns", call. = FALSE)
  }
}

#' @noRd
is_matrix <- function(x) (!is.character(x) && is.matrix(x)) || is(x, "sparseMatrix")

#' @noRd
validate_matrix <- function(x) {
  if (!is_matrix(x)) {
    name <- deparse(substitute(x))
    stop("`", name, "` must be a matrix or sparseMatrix", call. = FALSE)
  }
}
