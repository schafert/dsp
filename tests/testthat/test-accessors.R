# Fits used across the tests below. Kept small so the file runs quickly.

fit_cp <- suppressWarnings(suppressMessages(
  dsp_fit(y = {set.seed(200); c(rep(0, 40), rep(8, 40)) + rnorm(80)},
          model_spec = dsp_spec(family = "gaussian", model = "changepoint",
                                D = 1, useAnom = TRUE, obsSV = "SV"),
          nsave = 100, nburn = 100, verbose = FALSE)
))

sim_reg <- {set.seed(1); simRegression(nT = 40, p = 3, p_0 = 1)}

fit_reg <- suppressMessages(
  dsp_fit(y = sim_reg$y,
          model_spec = dsp_spec(family = "gaussian", model = "regression",
                                D = 1, X = sim_reg$X),
          nsave = 100, nburn = 100, verbose = FALSE)
)

fit_nb <- suppressMessages(
  dsp_fit(y = {set.seed(3); rpois(60, 5)},
          model_spec = dsp_spec(family = "negbinomial", model = "smoothing", D = 1),
          nsave = 100, nburn = 100, verbose = FALSE)
)


test_that("par_names lists the stored draws and nothing else", {

  nms <- par_names(fit_cp)

  expect_true(all(c("mu", "ypred", "omega", "gamma") %in% nms))
  expect_false(any(c("DIC", "p_d") %in% nms))

  # every name is accepted by the other accessors
  expect_equal(colnames(coda::as.mcmc(fit_cp, pars = "gamma")), "gamma")
  expect_error(par_names(list()), "class 'dsp'")
})


test_that("fitted returns pointwise summaries named as by quantile()", {

  f <- fitted(fit_cp)

  expect_equal(dim(f), c(80L, 6L))
  expect_equal(colnames(f), c("mean", "2.5%", "25%", "50%", "75%", "97.5%"))
  expect_equal(f[, "mean"], colMeans(fit_cp$mcmc_output$mu))

  expect_equal(colnames(fitted(fit_cp, probs = c(0.1, 0.9))), c("mean", "10%", "90%"))

  # a regression coefficient array keeps the predictor dimension
  b <- fitted(fit_reg, "beta")
  expect_equal(dim(b), c(40L, 3L, 6L))
  expect_equal(b[, 2, "mean"], colMeans(fit_reg$mcmc_output$beta[, , 2]))

  # btf_nb() returns r as a one-dimensional array; it is a scalar parameter
  r <- fitted(fit_nb, "r")
  expect_named(r, c("mean", "2.5%", "25%", "50%", "75%", "97.5%"))

  expect_error(fitted(fit_cp, c("mu", "ypred")), "single parameter")
  expect_error(fitted(fit_cp, "nope"), "must name one of")
})


test_that("DIC returns both estimators in a data frame", {

  d <- DIC(fit_cp)

  expect_s3_class(d, "data.frame")
  expect_equal(d$estimator, c("p_D", "p_V"))
  expect_equal(d$DIC, fit_cp$DIC$DIC)
  expect_equal(d$p_eff, fit_cp$DIC$p_d)

  no_dic <- fit_cp
  no_dic$DIC <- list()
  expect_error(DIC(no_dic), "computeDIC")
})


test_that("get_spec and nobs return the specification and series length", {

  expect_s3_class(get_spec(fit_cp), "dsp_spec")
  expect_identical(get_spec(fit_cp), fit_cp$model_spec)

  expect_equal(nobs(fit_cp), 80L)
  expect_equal(nobs(fit_reg), 40L)
})


test_that("as.mcmc is available without attaching coda", {

  expect_true("as.mcmc" %in% getNamespaceExports("dsp"))
  expect_identical(dsp::as.mcmc, coda::as.mcmc)
})


test_that("predict aligns the changepoint probabilities with the series", {

  cp <- predict(fit_cp, cp_thres = 0.5, cp_prop = TRUE)
  D <- fit_cp$model_spec$arguments$D

  expect_length(cp$cp_prop, 80L)
  expect_true(all(is.na(cp$cp_prop[seq_len(D)])))
  expect_true(all(cp$cp_prop[cp$cp_t] >= 0.5))
  expect_identical(cp$cp_t, predict(fit_cp, cp_thres = 0.5))
})


test_that("summary returns a summary.dsp object with one line per parameter", {

  s <- summary(fit_cp)

  expect_s3_class(s, "summary.dsp")
  expect_named(s, c("model_spec", "mcpar", "nobs", "DIC", "scalars", "vectors"))
  expect_equal(s$nobs, 80L)
  expect_identical(s$DIC, DIC(fit_cp))

  # gamma, dhs_phi, dhs_mean and loglike are scalars; the rest are vectors
  expect_setequal(rownames(s$scalars), c("gamma", "dhs_phi", "dhs_mean", "loglike"))
  expect_setequal(rownames(s$vectors), setdiff(par_names(fit_cp), rownames(s$scalars)))
  expect_equal(colnames(s$scalars), c("mean", "sd", "2.5%", "97.5%", "ess"))
  expect_equal(s$vectors["mu", "mean_max"], max(colMeans(fit_cp$mcmc_output$mu)))

  # pointwise effective sample sizes are opt in
  expect_null(s$vectors$ess_min)
  s_ess <- summary(fit_cp, pars = "mu", ess = TRUE)
  expect_true(s_ess$vectors["mu", "ess_min"] <= s_ess$vectors["mu", "ess_median"])

  # a regression coefficient array is reported one line per predictor
  expect_true(all(c("beta[, 1]", "beta[, 2]", "beta[, 3]") %in% rownames(summary(fit_reg)$vectors)))

  expect_true("r" %in% rownames(summary(fit_nb)$scalars))

  expect_warning(summary(fit_cp, pars = c("mu", "nope")), "nope")
  expect_error(summary(fit_cp, pars = "nope"), "None of the requested")
  expect_error(summary(fit_cp, probs = 0.5), "two probabilities")
})


test_that("the summary prints on one screen", {

  out <- capture.output(print(summary(fit_cp)))

  expect_lt(length(out), 30)
  expect_true(any(grepl("thinning interval 5", out)))
  expect_true(any(grepl("^DIC:", out)))
  capture.output(expect_invisible(print(summary(fit_cp))))
})


test_that("ess_varying excludes columns with no posterior variation", {

  b <- cbind(rep(1, 50), {set.seed(4); rnorm(50)})
  e <- ess_varying(b)

  expect_true(is.na(e[1]))
  expect_equal(e[2], unname(coda::effectiveSize(b[, 2])))
})
