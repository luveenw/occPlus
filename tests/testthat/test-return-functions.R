# verify all return*() and computeOccupancyProbs() functions correctly
# reshape the flat Stan draws matrix into labelled arrays.
#
# This is the highest-value tier: bugs here silently permute species or
# covariates in every downstream plot and analysis.
#
# sampleresults fixture: 1000 draws, 10 species, 100 sites, 3 primers, d=2,
# ncov_psi=2, ncov_theta=2, ncov_ord=2.
# Fixtures: loaded by helper-fixtures.R

# ---------------------------------------------------------------------------
# returnOccupancyCovariates()
# ---------------------------------------------------------------------------

test_that("returnOccupancyCovariates returns an array", {
  result <- returnOccupancyCovariates(fitmodel)
  expect_true(is.array(result))
})

test_that("returnOccupancyCovariates dimensions are [niter, ncov_psi, S]", {
  result <- returnOccupancyCovariates(fitmodel)
  expect_equal(dim(result), c(1000L, 2L, 10L))
})

test_that("returnOccupancyCovariates all values are finite", {
  result <- returnOccupancyCovariates(fitmodel)
  expect_true(all(is.finite(result)))
})

test_that("returnOccupancyCovariates covariate dimnames match X_psi columns", {
  result <- returnOccupancyCovariates(fitmodel)
  expect_equal(dimnames(result)[[2]], colnames(fitmodel$X_psi))
})

test_that("returnOccupancyCovariates species dimnames match infos$speciesNames", {
  result <- returnOccupancyCovariates(fitmodel)
  expect_equal(
    as.character(dimnames(result)[[3]]),
    as.character(fitmodel$infos$speciesNames)
  )
})

# ---------------------------------------------------------------------------
# returnDetectionCovariates()
# ---------------------------------------------------------------------------

test_that("returnDetectionCovariates returns an array", {
  result <- returnDetectionCovariates(fitmodel)
  expect_true(is.array(result))
})

test_that("returnDetectionCovariates dimensions are [niter, ncov_theta, S]", {
  result <- returnDetectionCovariates(fitmodel)
  expect_equal(dim(result), c(1000L, 2L, 10L))
})

test_that("returnDetectionCovariates all values are finite", {
  result <- returnDetectionCovariates(fitmodel)
  expect_true(all(is.finite(result)))
})

test_that("returnDetectionCovariates covariate dimnames match X_theta columns", {
  result <- returnDetectionCovariates(fitmodel)
  expect_equal(dimnames(result)[[2]], colnames(fitmodel$X_theta))
})

# ---------------------------------------------------------------------------
# returnOrdinationCovariatesOutput()
# ---------------------------------------------------------------------------

test_that("returnOrdinationCovariatesOutput returns an array", {
  result <- returnOrdinationCovariatesOutput(fitmodel)
  expect_true(is.array(result))
})

test_that("returnOrdinationCovariatesOutput dimensions are [niter, ncov_ord, d]", {
  result <- returnOrdinationCovariatesOutput(fitmodel)
  expect_equal(dim(result), c(1000L, fitmodel$infos$ncov_ord, fitmodel$infos$d))
})

test_that("returnOrdinationCovariatesOutput all values are finite", {
  result <- returnOrdinationCovariatesOutput(fitmodel)
  expect_true(all(is.finite(result)))
})

# ---------------------------------------------------------------------------
# returnOccupancyRates()
# ---------------------------------------------------------------------------

test_that("returnOccupancyRates returns an array", {
  result <- returnOccupancyRates(fitmodel)
  expect_true(is.array(result))
})

test_that("returnOccupancyRates dimensions are [niter, n_sites, S]", {
  result <- returnOccupancyRates(fitmodel)
  n_sites <- length(fitmodel$infos$siteNames)
  expect_equal(dim(result), c(1000L, n_sites, fitmodel$infos$S))
})

test_that("returnOccupancyRates all values are in [0, 1]", {
  result <- returnOccupancyRates(fitmodel)
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

# ---------------------------------------------------------------------------
# computeOccupancyProbs()
# ---------------------------------------------------------------------------

test_that("computeOccupancyProbs returns an array", {
  result <- computeOccupancyProbs(fitmodel)
  expect_true(is.array(result))
})

test_that("computeOccupancyProbs dimensions are [niter, n_sites, S]", {
  result <- computeOccupancyProbs(fitmodel)
  n_sites <- length(fitmodel$infos$siteNames)
  expect_equal(dim(result), c(1000L, n_sites, fitmodel$infos$S))
})

test_that("computeOccupancyProbs all values are in [0, 1]", {
  result <- computeOccupancyProbs(fitmodel)
  expect_true(all(result >= 0 & result <= 1, na.rm = TRUE))
})

test_that("computeOccupancyProbs site dimnames match infos$siteNames", {
  result <- computeOccupancyProbs(fitmodel)
  expect_equal(
    as.character(dimnames(result)[[2]]),
    as.character(fitmodel$infos$siteNames)
  )
})

test_that("computeOccupancyProbs species dimnames match infos$speciesNames", {
  result <- computeOccupancyProbs(fitmodel)
  expect_equal(
    as.character(dimnames(result)[[3]]),
    as.character(fitmodel$infos$speciesNames)
  )
})

# ---------------------------------------------------------------------------
# Consistency between returnOccupancyRates and computeOccupancyProbs
# ---------------------------------------------------------------------------

test_that("returnOccupancyRates and computeOccupancyProbs have same dimensions", {
  r1 <- returnOccupancyRates(fitmodel)
  r2 <- computeOccupancyProbs(fitmodel)
  expect_equal(dim(r1), dim(r2))
})

test_that("returnOccupancyRates and computeOccupancyProbs values agree within tolerance", {
  # Both compute logistic(X_psi * beta + U * LL) but via different code paths;
  # small floating-point differences are expected.
  r1 <- returnOccupancyRates(fitmodel)
  r2 <- computeOccupancyProbs(fitmodel)
  expect_equal(as.numeric(r1), as.numeric(r2), tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# Edge cases: idx_species subsetting
# ---------------------------------------------------------------------------

test_that("returnOccupancyCovariates values are consistent with direct draws extraction", {
  result <- returnOccupancyCovariates(fitmodel)
  # The first covariate, first species should match the raw draw for beta_psi[1,1]
  cn <- colnames(fitmodel$matrix_of_draws)
  raw_col <- fitmodel$matrix_of_draws[, cn == "beta_psi[1,1]"]
  expect_equal(as.numeric(result[, 1, 1]), as.numeric(raw_col), tolerance = 1e-10)
})

test_that("returnDetectionCovariates values are consistent with direct draws extraction", {
  result <- returnDetectionCovariates(fitmodel)
  cn <- colnames(fitmodel$matrix_of_draws)
  raw_col <- fitmodel$matrix_of_draws[, cn == "beta_theta[1,1]"]
  expect_equal(as.numeric(result[, 1, 1]), as.numeric(raw_col), tolerance = 1e-10)
})
