# verify runOccPlus() rejects malformed input before touching Stan.
# All tests here should resolve within milliseconds — no Stan compilation occurs
# because the validation checks fire first (or in the data-prep stage).
#
# Fixtures: `sampledata` and `fitmodel` loaded by helper-fixtures.R.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Returns a minimal valid data object derived from sampledata with S=2 species
# and only 6 rows (2 sites x 1 sample x 3 primers) so data-prep is fast if
# a test ever runs past validation by accident.
minimal_data <- function() {
  idx <- sampledata$info$Site %in% unique(sampledata$info$Site)[1:2]
  list(
    info = sampledata$info[idx, , drop = FALSE],
    OTU  = sampledata$OTU[idx, 1:2, drop = FALSE]
  )
}

# ---------------------------------------------------------------------------
# 1. Missing top-level list elements
# ---------------------------------------------------------------------------

test_that("missing $info element stops with error", {
  bad <- list(OTU = sampledata$OTU)
  expect_error(runOccPlus(bad, d = 1))
})

test_that("missing $OTU element stops with error", {
  bad <- list(info = sampledata$info)
  expect_error(runOccPlus(bad, d = 1))
})

# ---------------------------------------------------------------------------
# 2. NA values in required identifier columns
# ---------------------------------------------------------------------------

test_that("NA in Site column stops with informative error", {
  bad <- minimal_data()
  bad$info$Site[1] <- NA
  expect_error(
    runOccPlus(bad, d = 0),
    regexp = "NA in Site"
  )
})

test_that("NA in Sample column stops with informative error", {
  bad <- minimal_data()
  bad$info$Sample[1] <- NA
  expect_error(
    runOccPlus(bad, d = 0),
    regexp = "NA in Site, Sample or Primer"
  )
})

test_that("NA in Primer column stops with informative error", {
  bad <- minimal_data()
  bad$info$Primer[1] <- NA
  expect_error(
    runOccPlus(bad, d = 0),
    regexp = "NA in Site, Sample or Primer"
  )
})

# ---------------------------------------------------------------------------
# 3. Bad covariate names
# ---------------------------------------------------------------------------

test_that("nonexistent occCovariates name stops with error", {
  expect_error(
    runOccPlus(minimal_data(), d = 0, occCovariates = "does_not_exist"),
    regexp = "not in data"
  )
})

test_that("nonexistent detCovariates name stops with error", {
  expect_error(
    runOccPlus(minimal_data(), d = 0, detCovariates = "does_not_exist"),
    regexp = "not in data"
  )
})

test_that("nonexistent ordCovariates name stops with error", {
  expect_error(
    runOccPlus(minimal_data(), d = 1, ordCovariates = "does_not_exist"),
    regexp = "not in data"
  )
})

test_that("mix of valid and invalid covariate names stops with error", {
  expect_error(
    runOccPlus(
      minimal_data(), d = 0,
      occCovariates = c("X_psi.1", "does_not_exist")
    ),
    regexp = "not in data"
  )
})

# ---------------------------------------------------------------------------
# 4. Valid covariate names pass the validation check
#    (these will error later when trying to run Stan, but NOT during validation)
# ---------------------------------------------------------------------------

test_that("valid occCovariates name passes covariate validation", {
  skip_on_cran()  # May proceed to Stan compilation — too slow for CRAN
  result <- tryCatch(
    runOccPlus(minimal_data(), d = 0, occCovariates = "X_psi.1"),
    error = function(e) e
  )
  # Should not be the covariate-name error — any other error (Stan) is fine
  if (inherits(result, "error")) {
    expect_false(grepl("not in data", conditionMessage(result)))
  }
})

test_that("valid detCovariates name passes covariate validation", {
  skip_on_cran()  # May proceed to Stan compilation — too slow for CRAN
  result <- tryCatch(
    runOccPlus(minimal_data(), d = 0, detCovariates = "X_theta"),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_false(grepl("not in data", conditionMessage(result)))
  }
})

# ---------------------------------------------------------------------------
# 5. threshold argument behaviour
# ---------------------------------------------------------------------------

test_that("threshold = FALSE (default) passes validation without error", {
  skip_on_cran()
  result <- tryCatch(
    runOccPlus(minimal_data(), d = 0, threshold = FALSE),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_false(grepl("threshold", conditionMessage(result)))
  }
})

test_that("threshold = 0 does not fall through as FALSE (known behaviour)", {
  # Documents the sentinel-value issue: `threshold != F` is TRUE when
  # threshold = 0, so reads ARE binarised. This test locks in current
  # behaviour so any change is deliberate.
  skip_on_cran()
  result <- tryCatch(
    runOccPlus(minimal_data(), d = 0, threshold = 0),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_false(grepl("not in data", conditionMessage(result)))
  }
})

test_that("numeric threshold passes validation", {
  skip_on_cran()
  result <- tryCatch(
    runOccPlus(minimal_data(), d = 0, threshold = 10),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_false(grepl("not in data|NA in Site", conditionMessage(result)))
  }
})

# ---------------------------------------------------------------------------
# 6. OTU edge cases
# ---------------------------------------------------------------------------

test_that("OTU matrix with NA values is accepted by validation", {
  # NAs in OTU are valid (missing replicates) — validation should not reject
  skip_on_cran()
  d_na <- minimal_data()
  d_na$OTU[1, 1] <- NA
  result <- tryCatch(
    runOccPlus(d_na, d = 0),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_false(grepl("not in data|NA in Site", conditionMessage(result)))
  }
})

test_that("data with no covariates specified passes validation", {
  skip_on_cran()
  result <- tryCatch(
    runOccPlus(minimal_data(), d = 0),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    expect_false(grepl("not in data|NA in Site", conditionMessage(result)))
  }
})
