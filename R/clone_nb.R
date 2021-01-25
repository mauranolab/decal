#' @noRd
is_integer <- function(x) is.integer(x) || (is.numeric(x) && all.equal(x, round(x)))
#' @noRd
is_index <- function(x) is.integer(x) || (is.numeric(x) && all(x > 0) && all.equal(x, round(x)))

#' @noRd
as_index <- function(x, ref = NULL) {
  if (is_index(x) && (is.null(ref) || !is_integer(ref))) {
    x
  } else if (typeof(x) == typeof(ref)) {
    match(x, ref)
  } else {
    stop("`x` must either be an integer or match `ref` type")
    NULL
  }
}

#' @noRd
col_div <- function(mtx, vec) {
  if (length(vec) != ncol(mtx)) {
    stop("Incompatible dimensions at `col_div`. `vec` must have the same length as `mtx` columns.")
  }
  if (inherits(x = mtx, what = 'dgCMatrix')) {
    mtx@x <- mtx@x / vec[rep(seq_len(mtx@Dim[2]), diff(mtx@p))]
  } else if (inherits(x = mtx, what = 'dgTMatrix')) {
    mtx@x <- mtx@x / vec[mtx@j + 1]
  } else {
    mtx <- t(t(mtx) / vec)
  }
  return(mtx)
}

