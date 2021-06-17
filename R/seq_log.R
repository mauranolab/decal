#' Sequence generation in exponential scale
#'
#' Generates a numeric sequence with a specific length spaced equally in an
#' logarithm scale.
#'
#' @param from,to starting (minimal) and end (maximal) values fo the sequence
#' @param length_out desired length of the sequence
#' @param base logarithm scale base for which the data is spaced
#'
#' @export
seq_log <- function(from, to, length_out, base = 10L) {
  validate_numeric_scalar(from)
  validate_numeric_scalar(to)
  validate_positive_integer_scalar(length_out)
  validate_numeric_scalar(base)
  if (from >= to) {
    stop("`from` must be smaller than `to`", call. = FALSE)
  }

  log_seq <- seq(
    from = log(from, base = base),
    to = log(to, base = base),
    length.out = length_out
  )
  return(base**log_seq)
}
