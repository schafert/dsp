# Small fits shared by the tests below.
sim <- {set.seed(1); simRegression(nT = 40, p = 2, p_0 = 0)}
fit_reg <- suppressMessages(
  dsp_fit(sim$y, dsp_spec(family = "gaussian", model = "regression", D = 1,
                          X = sim$X),
          nsave = 100, nburn = 100, verbose = FALSE, seed = 1)
)
y_na <- {set.seed(2); cumsum(rnorm(50))}
y_na[10:12] <- NA
fit_s <- suppressMessages(
  dsp_fit(y_na, dsp_spec(family = "gaussian", model = "smoothing", D = 1),
          nsave = 100, nburn = 100, verbose = FALSE, seed = 1)
)

test_that("coef returns posterior mean coefficients named by predictor", {
  b <- coef(fit_reg)
  expect_equal(dim(b), c(40L, 2L))
  expect_equal(colnames(b), c("x1", "x2"))
  draws <- as.matrix(as.mcmc(fit_reg, pars = "beta"))
  expect_equal(unname(b[, 1]), unname(colMeans(draws)[1:40]))

  X_named <- sim$X
  colnames(X_named) <- c("market", "size")
  fit_named <- suppressMessages(
    dsp_fit(sim$y, dsp_spec(family = "gaussian", model = "regression", D = 1,
                            X = X_named),
            nsave = 20, nburn = 20, verbose = FALSE, seed = 1)
  )
  expect_equal(colnames(coef(fit_named)), c("market", "size"))

  expect_error(coef(fit_s), "use fitted\\(\\)")
})

test_that("residuals subtract the posterior mean and keep NA where y is missing", {
  r <- residuals(fit_s)
  expect_length(r, 50L)
  expect_true(all(is.na(r[10:12])))
  expect_equal(r[-(10:12)], (y_na - fitted(fit_s, pars = "mu")[, "mean"])[-(10:12)])
})

test_that("logLik returns a logLik object without a parameter count", {
  ll <- logLik(fit_s)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), mean(logLik(fit_s, summary = "draws")))
  expect_true(is.na(attr(ll, "df")))
  expect_equal(attr(ll, "nobs"), nobs(fit_s))
  expect_true(is.na(AIC(fit_s)))
  expect_true(is.na(BIC(fit_s)))
  expect_length(logLik(fit_s, summary = "draws"), 100L)
})

test_that("deviance is the posterior mean deviance used by DIC", {
  dev <- deviance(fit_s)
  expect_equal(dev, -2 * mean(logLik(fit_s, summary = "draws")))
  dic <- DIC(fit_s)
  expect_equal(dic$DIC[dic$estimator == "p_D"],
               dev + dic$p_eff[dic$estimator == "p_D"])
})

test_that("credible_intervals returns pointwise intervals labelled by time", {
  ci <- credible_intervals(fit_reg, pars = "beta", level = 0.8)
  expect_s3_class(ci, "data.frame")
  expect_equal(names(ci), c("time", "predictor", "mean", "lower", "upper"))
  expect_equal(nrow(ci), 80L)
  expect_equal(ci$time, rep(1:40, 2))
  expect_equal(ci$predictor, rep(1:2, each = 40))

  draws <- as.matrix(as.mcmc(fit_reg, pars = "beta"))
  expect_equal(c(ci$lower[1], ci$upper[1]),
               unname(quantile(draws[, 1], c(0.1, 0.9))))
  expect_equal(ci$mean[1:40], unname(coef(fit_reg)[, 1]))

  dates <- seq(as.Date("2000-01-01"), by = "month", length.out = 40)
  ci_mu <- credible_intervals(fit_reg, times = dates)
  expect_equal(names(ci_mu), c("time", "mean", "lower", "upper"))
  expect_equal(ci_mu$time, dates)
  fit_q <- fitted(fit_reg, probs = c(0.025, 0.975))
  expect_equal(ci_mu$lower, unname(fit_q[, "2.5%"]))
  expect_equal(ci_mu$upper, unname(fit_q[, "97.5%"]))

  # pointwise intervals lie within the simultaneous bands
  cb <- credible_bands(fit_reg, pars = "beta", level = 0.8)
  expect_true(all(ci$lower >= cb$lower - 1e-12 & ci$upper <= cb$upper + 1e-12))

  expect_error(credible_intervals(fit_reg, level = 1.2), "between 0 and 1")
  expect_error(credible_intervals(fit_reg, pars = "nope"), "one of the parameters")
  expect_error(credible_intervals(fit_reg, times = 1:3), "one entry per time point")
})

test_that("credible_intervals summarizes a scalar parameter in one row", {
  scalars <- Filter(function(p) is.null(dim(dsp:::get_draws(fit_s)[[p]])),
                    par_names(fit_s))
  skip_if(length(scalars) == 0L)
  ci <- credible_intervals(fit_s, pars = scalars[1])
  expect_equal(names(ci), c("mean", "lower", "upper"))
  expect_equal(nrow(ci), 1L)
})

test_that("B-spline coefficients are labelled by coefficient, not time", {
  x <- seq(0, 1, length.out = 60)
  y <- {set.seed(3); sin(6 * x) + rnorm(60, sd = 0.2)}
  fit_bs <- suppressMessages(
    dsp_fit(y, dsp_spec(family = "gaussian", model = "bspline", times = x,
                        num_knots = 10, D = 2, evol_error = "DHS"),
            nsave = 100, nburn = 100, verbose = FALSE, seed = 1)
  )
  K <- length(coef(fit_bs))
  ci <- credible_intervals(fit_bs, pars = "beta")
  expect_equal(names(ci)[1], "coefficient")
  expect_equal(ci$coefficient, seq_len(K))
  cb <- credible_bands(fit_bs, pars = "beta")
  expect_equal(names(cb)[1], "coefficient")
  expect_output(print(cb), paste0("\\(", K, " coefficients\\)"))
  expect_equal(names(credible_intervals(fit_bs))[1], "time")
})
