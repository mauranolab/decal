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
