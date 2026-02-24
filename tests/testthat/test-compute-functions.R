# verify internal utility functions (logit, logistic) and derived
# computations (correlation matrix structure).
#
# These tests are entirely self-contained: no fixtures required.

# ---------------------------------------------------------------------------
# logistic()
# ---------------------------------------------------------------------------

test_that("logistic(0) == 0.5 (symmetry point)", {
  expect_equal(logistic(0), 0.5)
})

test_that("logistic(-Inf) == 0 (lower boundary)", {
  expect_equal(logistic(-Inf), 0)
})

test_that("logistic(Inf) == 1 (upper boundary)", {
  expect_equal(logistic(Inf), 1)
})

test_that("logistic is monotonically increasing", {
  x <- seq(-5, 5, by = 0.5)
  y <- logistic(x)
  expect_true(all(diff(y) > 0))
})

test_that("logistic is vectorised over numeric vectors", {
  result <- logistic(c(-1, 0, 1))
  expect_length(result, 3L)
  expect_true(all(is.numeric(result)))
})

test_that("logistic is vectorised over matrices", {
  m <- matrix(c(-2, -1, 0, 1), nrow = 2)
  result <- logistic(m)
  expect_equal(dim(result), c(2L, 2L))
  expect_true(all(result > 0 & result < 1))
})

# ---------------------------------------------------------------------------
# logit()
# ---------------------------------------------------------------------------

test_that("logit(0.5) == 0 (symmetry point)", {
  expect_equal(logit(0.5), 0)
})

test_that("logit(0) == -Inf (lower boundary)", {
  expect_equal(logit(0), -Inf)
})

test_that("logit(1) == Inf (upper boundary)", {
  expect_equal(logit(1), Inf)
})

test_that("logit is vectorised", {
  result <- logit(c(0.1, 0.5, 0.9))
  expect_length(result, 3L)
})

# ---------------------------------------------------------------------------
# Round-trip: logistic(logit(x)) == x and logit(logistic(x)) == x
# ---------------------------------------------------------------------------

test_that("logistic(logit(p)) == p for p in (0, 1)", {
  p_vals <- c(0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
  expect_equal(logistic(logit(p_vals)), p_vals, tolerance = 1e-12)
})

test_that("logit(logistic(x)) == x for x in reals", {
  x_vals <- c(-5, -2, -0.5, 0, 0.5, 2, 5)
  expect_equal(logit(logistic(x_vals)), x_vals, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# generateCorrelationMatrixOutput() structural properties
# ---------------------------------------------------------------------------

test_that("correlation matrix has unit diagonal (all 1s at posterior median)", {
  Lambda <- generateCorrelationMatrixOutput(fitmodel)
  # Take the posterior median across iterations (dim 1)
  Lambda_med <- apply(Lambda, c(2, 3), median)
  expect_equal(diag(Lambda_med), rep(1, nrow(Lambda_med)), tolerance = 1e-10)
})

test_that("correlation matrix is symmetric at every posterior iteration", {
  Lambda <- generateCorrelationMatrixOutput(fitmodel)
  # Check a random sample of iterations
  set.seed(42)
  iter_check <- sample(dim(Lambda)[1], 10)
  for (iter in iter_check) {
    m <- Lambda[iter, , ]
    expect_equal(m, t(m), tolerance = 1e-10)
  }
})

test_that("correlation values are in [-1, 1]", {
  Lambda <- generateCorrelationMatrixOutput(fitmodel)
  expect_true(all(Lambda >= -1 - 1e-10 & Lambda <= 1 + 1e-10))
})

test_that("correlation matrix idx_species subset returns correct dimension", {
  Lambda_full <- generateCorrelationMatrixOutput(fitmodel)
  Lambda_sub  <- generateCorrelationMatrixOutput(fitmodel, idx_species = 1:3)
  # Full: niter x S x S; subset: niter x 3 x 3
  expect_equal(dim(Lambda_full)[1], dim(Lambda_sub)[1])  # same niter
  expect_equal(dim(Lambda_sub)[2], 3L)
  expect_equal(dim(Lambda_sub)[3], 3L)
})
