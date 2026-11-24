test_that("uni.slice errors instead of looping when the log density at x0 is -Inf", {
  # A collapsed scale parameter: the log density underflows to -Inf at x0.
  # Before the fix the stepping-out loop for R never ended.
  nT <- 229
  g <- function(x) -(nT + 2) * log(x) - 0.5 * 0.5 / x^2 - log(1 + (5 / x)^2)
  expect_identical(g(1e-160), -Inf)
  elapsed <- system.time(
    expect_error(uni.slice(1e-160, g, lower = 0, upper = Inf), "not finite")
  )[["elapsed"]]
  expect_lt(elapsed, 5)
})

test_that("uni.slice errors when the log density at x0 is NaN or +Inf", {
  expect_error(uni.slice(1, function(x) NaN), "not finite")
  expect_error(uni.slice(1, function(x) Inf), "not finite")
})

test_that("uni.slice caps the stepping-out and shrinkage phases", {
  # Flat log density on an unbounded support: stepping out never stops.
  expect_error(uni.slice(0, function(x) 0, max_iter = 100), "stepping out")
  # Slice of measure zero: proposals are rejected until the interval collapses
  # onto x0, which takes far more than 10 shrinkage steps.
  g_spike <- function(x) if (x == 0.5) 0 else -1e10
  expect_error(uni.slice(0.5, g_spike, lower = 0, upper = 1, max_iter = 10),
               "shrinkage")
})

test_that("uni.slice still samples the target distribution", {
  set.seed(1)
  x <- numeric(5000); x[1] <- 0
  for (i in 2:5000) x[i] <- uni.slice(x[i - 1], function(z) dnorm(z, log = TRUE))
  expect_lt(abs(mean(x)), 0.1)
  expect_lt(abs(sd(x) - 1), 0.1)
})
