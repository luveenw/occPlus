# Internal utility functions shared across output.R, output_mcmc.R, and mcmcfun.R.
# Single canonical definitions — do not redefine these in other files.

#' @keywords internal
logit <- function(x) {
  log(x / (1 - x))
}

#' @keywords internal
logistic <- function(x) {
  1 / (1 + exp(-x))
}
