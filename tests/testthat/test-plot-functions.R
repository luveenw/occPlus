# verify all plot*() functions return valid ggplot objects without error.
# Tests check class, layer count, and forced rendering via ggplot_build().
# They do NOT check visual appearance (no vdiffr snapshots).
#
# Fixtures: loaded by helper-fixtures.R

# ---------------------------------------------------------------------------
# Helper: force full ggplot render without opening a graphics device
# ---------------------------------------------------------------------------
renders_without_error <- function(p) {
  !inherits(tryCatch(ggplot2::ggplot_build(p), error = function(e) e), "error")
}

# ---------------------------------------------------------------------------
# Rate plots (no covariate argument required)
# ---------------------------------------------------------------------------

test_that("plotOccupancyRates returns a ggplot for a species subset", {
  p <- plotOccupancyRates(fitmodel, idx_species = 1:5)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotOccupancyRates returns a ggplot for all species", {
  p <- plotOccupancyRates(fitmodel, idx_species = 1:10)
  expect_s3_class(p, "ggplot")
  expect_gte(length(p$layers), 1L)
})

test_that("plotCollectionRates returns a ggplot", {
  p <- plotCollectionRates(fitmodel, idx_species = 1:5)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotDetectionRates returns a ggplot", {
  p <- plotDetectionRates(fitmodel, idx_species = 1:5)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotStage1FPRates returns a ggplot", {
  p <- plotStage1FPRates(fitmodel, idx_species = 1:5)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotStage2FPRates returns a ggplot", {
  p <- plotStage2FPRates(fitmodel, idx_species = 1:5)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

# ---------------------------------------------------------------------------
# Read intensity plot
# ---------------------------------------------------------------------------

test_that("plotReadIntensity returns a ggplot", {
  p <- plotReadIntensity(fitmodel)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotReadIntensity has at least one layer", {
  p <- plotReadIntensity(fitmodel)
  expect_gte(length(p$layers), 1L)
})

# ---------------------------------------------------------------------------
# Correlation matrix plots
# ---------------------------------------------------------------------------

test_that("plotCorrelationMatrix returns a ggplot", {
  # ggcorrplot uses the deprecated aes_string() — suppress that upstream warning
  p <- suppressWarnings(plotCorrelationMatrix(fitmodel))
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotSigElementsCorMatrix returns a ggplot", {
  p <- suppressWarnings(plotSigElementsCorMatrix(fitmodel))
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

# ---------------------------------------------------------------------------
# Covariate coefficient plots
# ---------------------------------------------------------------------------

test_that("plotOccupancyCovariates returns a ggplot for first covariate", {
  p <- plotOccupancyCovariates(fitmodel, covName = "X_psi.1", idx_species = 1:5)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotOccupancyCovariates returns a ggplot for second covariate", {
  p <- plotOccupancyCovariates(fitmodel, covName = "X_psi.2", idx_species = 1:5)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotOccupancyCovariates with NULL covName stops with informative error", {
  expect_error(
    plotOccupancyCovariates(fitmodel, covName = NULL),
    regexp = "No name provided"
  )
})

test_that("plotOccupancyCovariates with bad covName stops with informative error", {
  expect_error(
    plotOccupancyCovariates(fitmodel, covName = "nonexistent_covariate"),
    regexp = "Covariate name not found"
  )
})

test_that("plotOrdinationCovariates returns a ggplot for first ordination covariate", {
  p <- plotOrdinationCovariates(fitmodel, covName = "X_ord.1")
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotOrdinationCovariates returns a ggplot for second ordination covariate", {
  p <- plotOrdinationCovariates(fitmodel, covName = "X_ord.2")
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotDetectionCovariates returns a ggplot for X_theta covariate", {
  # Valid covName is "X_theta" (column 2 of X_theta; column 1 is the intercept)
  p <- plotDetectionCovariates(fitmodel, covName = "X_theta", idx_species = 1:5)
  expect_s3_class(p, "ggplot")
  expect_true(renders_without_error(p))
})

test_that("plotDetectionCovariates with bad covName stops with informative error", {
  expect_error(
    plotDetectionCovariates(fitmodel, covName = "nonexistent_covariate"),
    regexp = "Covariate name not found"
  )
})

# ---------------------------------------------------------------------------
# All plots have at least one ggplot layer
# ---------------------------------------------------------------------------

test_that("all plot functions produce objects with at least one layer", {
  plots <- suppressWarnings(list(
    plotOccupancyRates(fitmodel, idx_species = 1:5),
    plotCollectionRates(fitmodel, idx_species = 1:5),
    plotDetectionRates(fitmodel, idx_species = 1:5),
    plotStage1FPRates(fitmodel, idx_species = 1:5),
    plotStage2FPRates(fitmodel, idx_species = 1:5),
    plotReadIntensity(fitmodel),
    plotCorrelationMatrix(fitmodel),   # upstream ggcorrplot aes_string() warning
    plotOccupancyCovariates(fitmodel, covName = "X_psi.1", idx_species = 1:5),
    plotOrdinationCovariates(fitmodel, covName = "X_ord.1"),
    plotDetectionCovariates(fitmodel, covName = "X_theta", idx_species = 1:5)
  ))
  for (p in plots) {
    expect_gte(length(p$layers), 1L)
  }
})
