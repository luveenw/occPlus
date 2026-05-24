## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Package size

The installed package is approximately 73 MB, primarily due to bundled Stan
model files (~53 MB in `inst/stan/`) and a pre-fitted example model object
(~15 MB in `data/`). The Stan files are required at runtime by the
`runOccPlus()` backend, and the example model object supports the package
vignettes without requiring users to run long MCMC chains during installation.

## Vignettes

Both vignettes use `eval = FALSE` for code chunks that fit the Bayesian model,
as MCMC sampling exceeds CRAN's check time limits. All code in these chunks is
syntactically valid and tested outside of `R CMD check`.

## Downstream dependencies

There are currently no downstream dependencies on this package.
