# Fits used across the tests below. Kept small so the file runs quickly.

# y is simulated before the call: inside it, set.seed() would run after
# dsp_fit() sets its seed and replace the sampler's stream
y_cp <- {set.seed(200); c(rep(0, 40), rep(8, 40)) + rnorm(80)}
spec_cp <- dsp_spec(family = "gaussian", model = "changepoint",
                    D = 1, useAnom = TRUE, obsSV = "SV")
fit_cp <- suppressWarnings(suppressMessages(
  dsp_fit(y = y_cp, model_spec = spec_cp,
          nsave = 100, nburn = 100, verbose = FALSE, seed = 1)
))

sim_reg <- {set.seed(1); simRegression(nT = 40, p = 3, p_0 = 1)}

fit_reg <- suppressMessages(
  dsp_fit(y = sim_reg$y,
          model_spec = dsp_spec(family = "gaussian", model = "regression",
                                D = 1, X = sim_reg$X),
          nsave = 100, nburn = 100, verbose = FALSE, seed = 1)
)

fit_nb <- suppressMessages(
  dsp_fit(y = {set.seed(3); rpois(60, 5)},
          model_spec = dsp_spec(family = "negbinomial", model = "smoothing", D = 1),
          nsave = 100, nburn = 100, verbose = FALSE, seed = 1)
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

  # gamma, dhs_phi, dhs_mean, loglike and the one initial evolution variance
  # (D = 1) are scalars; the rest are vectors
  expect_setequal(rownames(s$scalars),
                  c("evol_sigma_t2[1]", "gamma", "dhs_phi", "dhs_mean", "loglike"))
  expect_setequal(rownames(s$vectors),
                  c(setdiff(par_names(fit_cp), c(rownames(s$scalars), "evol_sigma_t2")),
                    "evol_sigma_t2[2:80]"))
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


test_that("credible_bands returns joint bands and SimBaS", {

  b <- credible_bands(fit_cp, pars = "omega")
  pw <- fitted(fit_cp, pars = "omega", probs = c(0.025, 0.975))

  expect_equal(colnames(b), c("mean", "lower", "upper", "simBaS"))
  expect_equal(nrow(b), nrow(pw))
  expect_equal(unname(b[, "mean"]), unname(pw[, "mean"]))

  # joint bands are at least as wide as the pointwise intervals
  expect_true(all(b[, "lower"] <= pw[, "2.5%"] + 1e-8))
  expect_true(all(b[, "upper"] >= pw[, "97.5%"] - 1e-8))

  expect_true(all(b[, "simBaS"] >= 0 & b[, "simBaS"] <= 1))

  # A run of 100 draws can mix poorly between placing the jump at t = 40 or 41,
  # so the evidence for the jump is checked on a longer run
  fit_long <- suppressWarnings(suppressMessages(
    dsp_fit(y = y_cp, model_spec = spec_cp,
            nsave = 1000, nburn = 1000, verbose = FALSE, seed = 1)
  ))
  bl <- credible_bands(fit_long, pars = "omega")

  # the jump at t = 41 is flagged. The model may instead place the jump at
  # t = 40 and explain y[40] as an anomaly, so accept either increment.
  expect_lt(min(bl[(40:41) - 1, "simBaS"]), 0.05)

  # the score is the smallest alpha at which the joint band excludes zero
  i <- which.min(bl[, "simBaS"])
  alpha <- bl[i, "simBaS"]
  excl <- credible_bands(fit_long, pars = "omega", level = 1 - alpha - 0.01)
  expect_true(excl[i, "lower"] > 0 || excl[i, "upper"] < 0)

  # one band per predictor
  bb <- credible_bands(fit_reg, pars = "beta")
  expect_equal(dim(bb), c(nobs(fit_reg), ncol(sim_reg$X), 4))
  expect_equal(dimnames(bb)[[3]], c("mean", "lower", "upper", "simBaS"))

  expect_error(credible_bands(fit_cp, pars = "gamma"), "scalar")
  expect_error(credible_bands(fit_cp, level = 95), "level")
  expect_error(credible_bands(fit_cp, pars = "nope"), "one of the parameters")
})


test_that("band helpers ignore columns with no posterior variation", {

  set.seed(1)
  m <- cbind(0, matrix(rnorm(400), 100, 4), 5)

  cb <- dsp:::credBands(m)
  expect_false(anyNA(cb))
  expect_equal(cb[1, ], c(0, 0))
  expect_equal(cb[6, ], c(5, 5))

  s <- dsp:::simBaS(m)
  expect_true(is.na(s[1]) && is.na(s[6]))
  expect_equal(s[2:5], dsp:::simBaS(m[, 2:5]))

  # agrees with the definition of Meyer et al. (2015)
  Efx <- colMeans(m[, 2:5]); SDfx <- apply(m[, 2:5], 2, sd)
  Maxfx <- apply(abs(sweep(sweep(m[, 2:5], 2, Efx), 2, SDfx, "/")), 1, max)
  expect_equal(s[2:5], sapply(abs(Efx)/SDfx, \(z) mean(z <= Maxfx)))
})


test_that("a constant observation variance is stored once per draw", {

  expect_null(dim(fit_reg$mcmc_output$obs_sigma_t2))
  expect_length(fit_reg$mcmc_output$obs_sigma_t2, 100)
  expect_gt(stats::sd(fit_reg$mcmc_output$obs_sigma_t2), 0)
  expect_length(fitted(fit_reg, pars = "obs_sigma_t2"), 6)

  # time-varying under SV
  expect_equal(dim(fit_cp$mcmc_output$obs_sigma_t2), c(100, 80))

  s <- summary(fit_reg)
  expect_true("obs_sigma_t2" %in% rownames(s$scalars))
  expect_false("obs_sigma_t2" %in% rownames(s$vectors))

  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(fit_reg, type = "obs_sigma_t2"))
})


test_that("summary reports the initial evolution variances separately", {

  s <- summary(fit_cp)
  # D = 1: one initial variance, a scalar
  expect_true("evol_sigma_t2[1]" %in% rownames(s$scalars))
  expect_true("evol_sigma_t2[2:80]" %in% rownames(s$vectors))
  expect_false("evol_sigma_t2" %in% c(rownames(s$scalars), rownames(s$vectors)))

  # regression: one pair of lines per predictor
  s <- summary(fit_reg)
  expect_true(all(paste0("evol_sigma_t2[1, ", 1:3, "]") %in% rownames(s$scalars)))
  expect_true(all(paste0("evol_sigma_t2[2:40, ", 1:3, "]") %in% rownames(s$vectors)))
})


test_that("the anomaly is aligned with y and zero for the first D", {

  expect_equal(dim(fit_cp$mcmc_output$zeta), c(100, 80))
  expect_equal(dim(fit_cp$mcmc_output$zeta_sigma_t2), c(100, 80))
  # D = 1: the anomaly is defined for t > D only
  expect_true(all(fit_cp$mcmc_output$zeta[, 1] == 0))
  expect_true(all(fit_cp$mcmc_output$zeta_sigma_t2[, 1] == 0))
  expect_true(all(fit_cp$mcmc_output$zeta_sigma_t2[, -1] > 0))
  expect_equal(dim(fit_cp$mcmc_output$omega), c(100, 79))

  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(fit_cp, type = "zeta_sigma_t2"))
})
