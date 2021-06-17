#' @importFrom Matrix t
#' @noRd
coldiv <- function(mtx, vec) {
  ## Validate input size
  if (length(vec) != ncol(mtx)) {
    stop(
      "Incompatible dimensions. `vec` must be the same size as `mtx` columns",
      call. = FALSE
    )
  }
  if (inherits(x = mtx, what = "dgCMatrix")) {
    mtx@x <- mtx@x / vec[rep(seq_len(mtx@Dim[2]), diff(mtx@p))]
  } else if (inherits(x = mtx, what = "dgTMatrix")) {
    mtx@x <- mtx@x / vec[mtx@j + 1]
  } else if (inherits(x = mtx, what = "Matrix")) {
    mtx <- Matrix::t(Matrix::t(mtx) / vec)
  } else {
    mtx <- t(t(mtx) / vec)
  }
  return(mtx)
}

#' @noRd
or <- function(x, default) {
  if (is.null(x)) {
    return(default)
  }
  return(x)
}

#' @noRd
get_index <- function(x, ref = NULL) {
  if (!is.null(ref) && typeof(x) == typeof(ref)) {
    return(match(x, ref))
  }
  if (is.integer(x) && all(x > 0L)) {
    return(x)
  }
  stop("`x` must be an integer index or match reference type")
}
