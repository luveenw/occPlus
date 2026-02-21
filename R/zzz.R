#' @useDynLib occPlus, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  # Set mc.cores (but limit to 2 for CRAN compliance)
  options(mc.cores = min(2L, parallel::detectCores()))
  # Cache compiled Stan models to disk to avoid recompilation on each session
  rstan::rstan_options(auto_write = TRUE)
}


