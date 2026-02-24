# specification tests for the custom MCMC backend (runMCMCOccPlus).
#
# ALL TESTS ARE SKIPPED because the MCMC backend depends on Rcpp functions
# compiled from src/functions.cpp, which is not present in this repository.
#
# Purpose: these tests serve as a machine-readable specification.  When
# src/functions.cpp is provided and the MCMC backend becomes operational,
# unskip these tests by removing the skip() calls.  The tests will then
# immediately verify whether the MCMC output is structurally compatible with
# the Stan output and whether all downstream output functions accept it.
#
# Fixtures: loaded by helper-fixtures.R

SKIP_MSG <- paste(
  "MCMC backend requires src/functions.cpp (Rcpp) which is",
  "not present in this repository"
)

# Minimal dataset for a fast MCMC run (2 sites, 1 sample/site, 2 primers,
# 2 species, nchain=1, nburn=10, niter=20).
mcmc_minimal_data <- local({
  idx <- sampledata$info$Site %in% unique(sampledata$info$Site)[1:2]
  list(
    info = sampledata$info[idx, , drop = FALSE],
    OTU  = sampledata$OTU[idx, 1:2, drop = FALSE]
  )
})

mcmc_params <- list(nchain = 1L, nburn = 10L, niter = 20L)

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("runMCMCOccPlus returns a list with required elements", {
  skip(SKIP_MSG)
  result <- runMCMCOccPlus(
    mcmc_minimal_data, d = 1L,
    MCMCparams = mcmc_params
  )
  expect_type(result, "list")
  expect_named(result, c("results_output", "infos", "X_ord", "X_theta", "X_psi"),
               ignore.order = TRUE)
})

test_that("runMCMCOccPlus results_output contains all expected arrays", {
  skip(SKIP_MSG)
  result <- runMCMCOccPlus(
    mcmc_minimal_data, d = 1L,
    MCMCparams = mcmc_params
  )
  ro <- result$results_output
  expected_names <- c(
    "beta_ord_output", "beta_psi_output", "beta_theta_output",
    "LL_output", "E_output", "U_output", "z_output",
    "p_output", "q_output", "theta0_output"
  )
  for (nm in expected_names) {
    expect_true(nm %in% names(ro),
                label = paste("results_output contains", nm))
  }
})

test_that("runMCMCOccPlus infos element has same fields as Stan infos", {
  skip(SKIP_MSG)
  result <- runMCMCOccPlus(
    mcmc_minimal_data, d = 1L,
    MCMCparams = mcmc_params
  )
  stan_info_names <- names(fitmodel$infos)
  mcmc_info_names <- names(result$infos)
  for (nm in stan_info_names) {
    expect_true(nm %in% mcmc_info_names,
                label = paste("MCMC infos contains field:", nm))
  }
})

# ---------------------------------------------------------------------------
# Array dimensions
# ---------------------------------------------------------------------------

test_that("beta_psi_output has 4 dimensions [ncov_psi, S, niter, nchain]", {
  skip(SKIP_MSG)
  result <- runMCMCOccPlus(
    mcmc_minimal_data, d = 1L,
    MCMCparams = mcmc_params
  )
  bp <- result$results_output$beta_psi_output
  expect_equal(length(dim(bp)), 4L)
  # dim[3] == niter, dim[4] == nchain
  expect_equal(dim(bp)[3], mcmc_params$niter)
  expect_equal(dim(bp)[4], mcmc_params$nchain)
})

test_that("z_output has 4 dimensions [n_sites, S, niter, nchain]", {
  skip(SKIP_MSG)
  result <- runMCMCOccPlus(
    mcmc_minimal_data, d = 1L,
    MCMCparams = mcmc_params
  )
  z <- result$results_output$z_output
  expect_equal(length(dim(z)), 4L)
  expect_equal(dim(z)[3], mcmc_params$niter)
})

# ---------------------------------------------------------------------------
# Compatibility with output functions
# ---------------------------------------------------------------------------

test_that("returnOrdinationCovariatesOutput_mcmc works with MCMC result", {
  skip(SKIP_MSG)
  result <- runMCMCOccPlus(
    mcmc_minimal_data, d = 1L,
    MCMCparams = mcmc_params
  )
  out <- returnOrdinationCovariatesOutput_mcmc(result)
  expect_true(is.array(out))
  expect_true(all(is.finite(out)))
})
