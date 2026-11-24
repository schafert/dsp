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

test_that("confint uses the generic's parm argument and accepts pars", {
  ci <- confint(fit_reg, parm = "beta", level = 0.8)
  expect_equal(colnames(ci), c("10%", "90%"))
  expect_equal(nrow(ci), 80L)
  expect_true(all(grepl("^beta", rownames(ci))))
  expect_identical(confint(fit_reg, pars = "beta", level = 0.8), ci)
  expect_identical(confint(fit_reg, "beta", level = 0.8), ci)
  expect_true(all(grepl("^mu", rownames(confint(fit_reg)))))

  draws <- as.matrix(as.mcmc(fit_reg, pars = "beta"))
  expect_equal(unname(ci[1, ]), unname(quantile(draws[, 1], c(0.1, 0.9))))

  expect_error(confint(fit_reg, parm = "beta", pars = "mu"), "not both")
  expect_error(confint(fit_reg, parameter = "beta"), "Unused")
  expect_error(confint(fit_reg, level = 1.2), "between zero and one")
})
