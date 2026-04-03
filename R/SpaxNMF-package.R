#' SpaxNMF: Regularized Non-negative Matrix Factorization for IFU Cubes
#'
#' Tools for reshaping IFU cubes, fitting smooth non-negative factorization
#' models, and building reproducible examples for package documentation and
#' comparison experiments.
#'
#' The package is designed around a simple workflow:
#' \enumerate{
#'   \item simulate or ingest a spectral cube;
#'   \item reshape it with [cube_to_matrix()];
#'   \item fit a low-rank unmixing model with [spax_nmf()];
#'   \item inspect component spectra and component-weight maps.
#' }
#'
#' @docType package
#' @name SpaxNMF
NULL
