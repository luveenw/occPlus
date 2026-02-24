# verify the data-preparation products stored in the pre-computed
# sampleresults fixture.  These tests do not call Stan — they inspect the
# design matrices, infos list, and posterior draws matrix that runOccPlus()
# already produced and saved.
#
# sampledata: 100 sites, 3 samples/site, 3 primers/sample, 10 species, d=2,
#             occCovariates = c("X_psi.1","X_psi.2"),
#             ordCovariates = c("X_ord.1","X_ord.2"),
#             detCovariates = c("X_theta")
# Fixtures: loaded by helper-fixtures.R

# ---------------------------------------------------------------------------
# infos metadata
# ---------------------------------------------------------------------------

test_that("infos$S equals number of species (10)", {
  expect_equal(fitmodel$infos$S, 10L)
})

test_that("infos$d equals number of latent factors (2)", {
  expect_equal(fitmodel$infos$d, 2L)
})

test_that("infos$ncov_psi equals number of occupancy covariate columns", {
  expect_equal(fitmodel$infos$ncov_psi, 2L)
})

test_that("infos$ncov_theta equals number of detection covariate columns", {
  expect_equal(fitmodel$infos$ncov_theta, 2L)
})

test_that("infos$ncov_ord equals number of ordination covariate columns", {
  expect_equal(fitmodel$infos$ncov_ord, 2L)
})

test_that("infos$siteNames has 100 entries (one per site)", {
  expect_length(fitmodel$infos$siteNames, 100L)
})

test_that("infos$primerNames has 3 entries", {
  expect_length(fitmodel$infos$primerNames, 3L)
})

test_that("infos$speciesNames has 10 entries", {
  expect_length(fitmodel$infos$speciesNames, 10L)
})

# ---------------------------------------------------------------------------
# Design matrix dimensions and column names
# ---------------------------------------------------------------------------

test_that("X_psi has correct dimensions: n_sites x ncov_psi", {
  expect_equal(dim(fitmodel$X_psi), c(100L, 2L))
})

test_that("X_psi column names match occupancy covariate names", {
  expect_equal(colnames(fitmodel$X_psi), c("X_psi.1", "X_psi.2"))
})

test_that("X_theta has correct dimensions: n_samples x ncov_theta", {
  # 100 sites * 3 samples/site = 300 samples; ncov_theta = 2 (intercept + X_theta)
  expect_equal(dim(fitmodel$X_theta), c(300L, 2L))
})

test_that("X_theta has column names", {
  expect_false(is.null(colnames(fitmodel$X_theta)))
  expect_length(colnames(fitmodel$X_theta), 2L)
})

test_that("X_ord has correct dimensions: n_sites x ncov_ord", {
  expect_equal(dim(fitmodel$X_ord), c(100L, 2L))
})

test_that("X_ord column names match ordination covariate names", {
  expect_equal(colnames(fitmodel$X_ord), c("X_ord.1", "X_ord.2"))
})

# ---------------------------------------------------------------------------
# Posterior draws matrix
# ---------------------------------------------------------------------------

test_that("matrix_of_draws has 1000 rows (posterior samples)", {
  expect_equal(nrow(fitmodel$matrix_of_draws), 1000L)
})

test_that("matrix_of_draws has no NA values", {
  expect_false(any(is.na(fitmodel$matrix_of_draws)))
})

test_that("matrix_of_draws has column names", {
  expect_false(is.null(colnames(fitmodel$matrix_of_draws)))
})

test_that("matrix_of_draws contains expected parameter prefixes", {
  cn <- colnames(fitmodel$matrix_of_draws)
  expect_true(any(startsWith(cn, "beta_psi[")))
  expect_true(any(startsWith(cn, "beta_theta[")))
  expect_true(any(startsWith(cn, "beta_ord[")))
  expect_true(any(startsWith(cn, "beta0_psi[")))
  expect_true(any(startsWith(cn, "U[")))
  expect_true(any(startsWith(cn, "LL[")))
  expect_true(any(startsWith(cn, "p[")))
  expect_true(any(startsWith(cn, "q[")))
  expect_true(any(startsWith(cn, "theta0[")))
  expect_true(any(cn == "mu1"))
  expect_true(any(cn == "sigma0"))
  expect_true(any(cn == "sigma1"))
  expect_true(any(cn == "pi0"))
  expect_true(any(cn == "lp__"))
})

test_that("matrix_of_draws beta_psi count matches ncov_psi * S", {
  cn <- colnames(fitmodel$matrix_of_draws)
  n_beta_psi <- sum(startsWith(cn, "beta_psi["))
  expect_equal(n_beta_psi, fitmodel$infos$ncov_psi * fitmodel$infos$S)
})

test_that("matrix_of_draws beta_theta count matches ncov_theta * S", {
  cn <- colnames(fitmodel$matrix_of_draws)
  n_beta_theta <- sum(startsWith(cn, "beta_theta["))
  expect_equal(n_beta_theta, fitmodel$infos$ncov_theta * fitmodel$infos$S)
})

test_that("matrix_of_draws U count matches n_sites * d", {
  cn <- colnames(fitmodel$matrix_of_draws)
  n_U <- sum(startsWith(cn, "U["))
  n_sites <- length(fitmodel$infos$siteNames)
  expect_equal(n_U, n_sites * fitmodel$infos$d)
})

# ---------------------------------------------------------------------------
# Stan fit object
# ---------------------------------------------------------------------------

test_that("vb_fit is a stanfit S4 object", {
  expect_s4_class(fitmodel$vb_fit, "stanfit")
})

test_that("vb_fit@model_name is not empty", {
  expect_true(nchar(fitmodel$vb_fit@model_name) > 0)
})
