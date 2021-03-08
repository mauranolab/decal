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

  log_seq <- seq(
    from = log(from, base = base),
    to = log(to, base = base),
    length.out = length_out
  )
  return(base**log_seq)
}

#' @noRd
get_ordered <- function(x, order = NULL) {
  if (!is.null(order)) {
    notin <- !x %in% order
    if (any(notin)) order <- c(order, sort(unique(x[notin])))
    return(order)
  }
  if (is.numeric(x) && x == round(x)) {
    rg <- range(x)
    return(seq(rg[1], rg[2]))
  }
  return(sort(unique(x)))
}

#' Builds a adjacency matrix of clone and cells
#'
#' Given the cell composition of each clone, it builds and returns a binary
#' adjacency sparse matrix where each clone corresponds to a column, each
#' cell to a row and a 1 value indicates that cell belongs to the clone.
#'
#' @param clones vector of clone ids
#' @param cells vector or list of cell ids with the same length as `clones`
#' @param cell_order ordered vector of cell ids indicating their row position.
#' When `NULL` it uses the sorted order of elements,
#' @return a sparse adjacency matrix indicating which cell belongs to which
#' clone.
#'
#' @importFrom Matrix sparseMatrix
#' @export
as_clone_matrix <- function(clones, cells, cell_order = NULL) {
  if (length(clones) != length(cells)) {
    stop("Both `clones` and `cells` must have the same length", call. = FALSE)
  }
  if (is.list(cells)) {
    clones <- rep(clones, sapply(cells, length))
    cells <- unlist(cells)
  }
  validate_index(clones)
  validate_index(cells)
  clone_order <- get_ordered(clones)
  cell_order <- get_ordered(cells, cell_order)

  clone_matrix <- sparseMatrix(
    i = match(cells, cell_order),
    j = match(clones, clone_order),
    x = 1L,
    dims = c(length(cell_order), length(clone_order))
  )

  if (is.character(clone_order)) colnames(clone_matrix) <- clone_order
  if (is.character(cell_order)) rownames(clone_matrix) <- cell_order

  return(clone_matrix)
}
