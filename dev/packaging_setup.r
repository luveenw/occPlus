# Install dependencies - this is already taken care of by the DESCRIPTION file, but you can run this to make sure everything is installed
install.packages(c("fansi", "downlit", "pkgdown", "devtools", "roxygen2", "usethis"))
# pak::pak(c("fansi", "downlit", "pkgdown", "devtools", "roxygen2", "usethis"))

# One-time calls to add Rcpp and Armadillo as dependencies for C++ compilation
usethis::use_package("Rcpp", type = "LinkingTo")
usethis::use_package("RcppArmadillo", type = "LinkingTo")

# one-time call to check documentation
# devtools::check(args = c("--as-cran"))

# install and build from GitHub (WIP)
# remotes::install_github(
#   "luveenw/occPlus@luveenw_test_plan_3", build_vignettes = FALSE)

# LOCAL ALTERNATIVE: if you downloaded the code, run this instead and comment
# out the line above this:
devtools::install(build_vignettes = FALSE)

# Build binary package
#
# NOTE: THIS GENERATES A PLATFORM-SPECIFIC PACKAGE. PACKAGE FILES GENERATED ON WINDOWS WON'T RUN ON MACOS, FOR EXAMPLE
#
# devtools::build(binary = TRUE)

library(occPlus)
getNamespaceExports("occPlus") # lists all exported functions and data
