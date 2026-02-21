.onLoad <- function(libname, pkgname) {
  # CRAN policy: packages must not set mc.cores > 2 without user consent.
  # We set it to the minimum of 2 and the detected core count.
  options(mc.cores = min(2L, parallel::detectCores()))
  # Cache compiled Stan models to disk to avoid recompilation on each session.
  rstan::rstan_options(auto_write = TRUE)
}


