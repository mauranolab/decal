#' Simulated scRNA-seq count matrix (1500 x 500)
#'
#' The matrix was generated using `sim_count()` for a total cell depth
#' following a log-normal distribution (\eqn{ln(D) ~ N(9, .5)}), a gene
#' expression rate ranging from 0.0001 to 0.01 in a logarithm scale
#' and overdispersion (`theta`) of 80.
#'
#' @format A sparse count matrix with 1500 rows (features) and 500 columns
#' (cells).
"sc_simulated"
