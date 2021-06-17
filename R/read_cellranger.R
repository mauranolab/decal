#' @noRd
check_file <- function(...) {
  file <- file.path(...)
  if (!file.exists(file)) {
    stop("Missing required file `", file, "`", call. = FALSE)
  }
  return(file)
}

#' @noRd
read_tsv_ <- function(file, ...) {
  read.delim(file,
    header = FALSE, sep = "\t", row.names = NULL,
    stringsAsFactors = FALSE, ...
  )
}

#' Load data from cellranger project
#'
#' Enables easy loading of 10x single-cell experiment sparse matrix data.
#' Currently compatible only with Cell Ranger >= v3.0
#'
#' @param data_dir Directory containing matrix.mtx, features.tsv, and
#' barcodes.tsv file produced by 10X run.
#' @return a sparse matrix containing the expression data.
#'
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
#' @export
read_cellranger <- function(data_dir) {
  ## Largely based on `seurat::Read10X`
  if (!is.character(data_dir) || length(data_dir) != 1) {
    stop("Invalid `data_dir` argument", call. = FALSE)
  }
  if (!dir.exists(data_dir)) {
    stop("Unable to find `", data_dir, "`", call. = FALSE)
  }
  feature_file <- check_file(data_dir, "features.tsv.gz")
  barcode_file <- check_file(data_dir, "barcodes.tsv.gz")
  matrix_file <- check_file(data_dir, "matrix.mtx.gz")

  sparse_dat <- readMM(file = matrix_file)
  colnames(sparse_dat) <- read_tsv_(barcode_file)[, 1]
  rownames(sparse_dat) <- read_tsv_(feature_file)[, 1]
  return(sparse_dat)
}
