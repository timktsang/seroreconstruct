#' @keywords internal
"_PACKAGE"

#' @useDynLib seroreconstruct, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom grDevices adjustcolor hcl.colors rgb
#' @importFrom graphics layout par plot abline axis box lines points legend polygon segments
#' @importFrom stats quantile density median acf
NULL

.onLoad <- function(libname, pkgname) {
  # CRAN policy: packages must not use more than 2 threads during checks.
  # Users can override via RCPP_PARALLEL_NUM_THREADS environment variable.
  RcppParallel::setThreadOptions(numThreads = 2L)
}
