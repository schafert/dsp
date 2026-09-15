# Fits used across the tests below. Kept small so the file runs quickly.

fit_smooth <- suppressMessages(
  dsp_fit(y = {set.seed(200); c(rep(0, 50), rep(10, 50)) + rnorm(100)},
          model_spec = dsp_spec(family = "gaussian", model = "smoothing", D = 1),
          nsave = 100, nburn = 100, verbose = FALSE)
)

sim_reg <- {set.seed(1); simRegression(nT = 40, p = 3, p_0 = 1)}

fit_reg <- suppressMessages(
  dsp_fit(y = sim_reg$y,
          model_spec = dsp_spec(family = "gaussian", model = "regression",
                                D = 1, X = sim_reg$X),
          nsave = 100, nburn = 100, verbose = FALSE)
)


test_that("get_flatnames indexes in column-major order", {

  # A scalar parameter keeps its bare name
  expect_equal(get_flatnames("dhs_phi"), "dhs_phi")
  expect_equal(get_flatnames("dhs_phi", integer(0)), "dhs_phi")

  # A vector parameter is indexed once
  expect_equal(get_flatnames("mu", 4),
               c("mu[1]", "mu[2]", "mu[3]", "mu[4]"))

  # A matrix parameter varies the first index fastest, as as.vector() does
  expect_equal(get_flatnames("beta", c(3, 2)),
               c("beta[1, 1]", "beta[2, 1]", "beta[3, 1]",
                 "beta[1, 2]", "beta[2, 2]", "beta[3, 2]"))

  expect_length(get_flatnames("beta", c(80, 3)), 240)

  # The names must agree with the flattening they label
  a <- array(1:24, c(2, 3, 4))
  flat <- matrix(as.vector(a), nrow = 2)
  nms <- get_flatnames("b", dim(a)[-1])
  expect_equal(nms[6], "b[3, 2]")
  expect_equal(as.numeric(flat[, 6]), as.numeric(a[, 3, 2]))

  expect_equal(get_flatnames("mu", 3, first_index = 0),
               c("mu[0]", "mu[1]", "mu[2]"))
})


test_that("as.mcmc.dsp returns an mcmc object with the sampler's iteration index", {

  draws <- coda::as.mcmc(fit_smooth)

  expect_s3_class(draws, "mcmc")
  expect_equal(nrow(draws), unname(fit_smooth$mcpar["nsave"]))

  # The sampler runs nburn + (nskip + 1) * nsave iterations and saves every
  # (nskip + 1) of them after burn-in, so the coda thinning interval is
  # nskip + 1 rather than nskip.
  nsave <- unname(fit_smooth$mcpar["nsave"])
  nburn <- unname(fit_smooth$mcpar["nburn"])
  nskip <- unname(fit_smooth$mcpar["nskip"])
  thin  <- nskip + 1

  expect_equal(attr(draws, "mcpar"),
               c(nburn + thin, nburn + nsave * thin, thin))
  expect_equal(coda::thin(draws), thin)
  expect_equal(stats::start(draws), nburn + thin)
  expect_equal(stats::end(draws), nburn + nsave * thin)
})


test_that("as.mcmc.dsp preserves the draws for vector and scalar parameters", {

  draws <- coda::as.mcmc(fit_smooth)

  expect_equal(as.numeric(draws[, "mu[7]"]),
               as.numeric(fit_smooth$mcmc_output$mu[, 7]))
  expect_equal(as.numeric(draws[, "mu[100]"]),
               as.numeric(fit_smooth$mcmc_output$mu[, 100]))
  expect_equal(as.numeric(draws[, "ypred[42]"]),
               as.numeric(fit_smooth$mcmc_output$ypred[, 42]))
  expect_equal(as.numeric(draws[, "loglike"]),
               as.numeric(fit_smooth$mcmc_output$loglike))
  expect_equal(as.numeric(draws[, "dhs_phi"]),
               as.numeric(fit_smooth$mcmc_output$dhs_phi))
})


test_that("as.mcmc.dsp preserves the draws for three-dimensional parameters", {

  draws <- coda::as.mcmc(fit_reg)

  beta_cols <- grep("^beta", colnames(draws), value = TRUE)
  expect_length(beta_cols, prod(dim(fit_reg$mcmc_output$beta)[-1]))
  expect_equal(head(beta_cols, 2), c("beta[1, 1]", "beta[2, 1]"))

  expect_equal(as.numeric(draws[, "beta[5, 2]"]),
               as.numeric(fit_reg$mcmc_output$beta[, 5, 2]))
  expect_equal(as.numeric(draws[, "evol_sigma_t2[9, 3]"]),
               as.numeric(fit_reg$mcmc_output$evol_sigma_t2[, 9, 3]))

  # Total width is the sum of the flattened parameter sizes
  expected_width <- sum(vapply(
    fit_reg$mcmc_output[!names(fit_reg$mcmc_output) %in% c("DIC", "p_d")],
    function(s) if (is.null(dim(s))) 1 else prod(dim(s)[-1]), numeric(1)))
  expect_equal(ncol(draws), expected_width)
})


test_that("as.mcmc.dsp drops posterior summaries rather than draws", {

  draws <- coda::as.mcmc(fit_smooth)

  expect_false(any(c("DIC", "p_d") %in% colnames(draws)))
  expect_false(any(grepl("^DIC|^p_d", colnames(draws))))

  # An element whose leading dimension is not nsave is not a monitored
  # quantity and must be filtered even if it is added later
  spiked <- fit_smooth
  spiked$mcmc_output$not_a_draw <- rnorm(7)
  expect_false(any(grepl("not_a_draw", colnames(coda::as.mcmc(spiked)))))
})


test_that("as.mcmc.dsp honours the pars argument", {

  draws <- coda::as.mcmc(fit_smooth, pars = c("dhs_phi", "dhs_mean"))

  expect_equal(colnames(draws), c("dhs_phi", "dhs_mean"))
  expect_equal(nrow(draws), unname(fit_smooth$mcpar["nsave"]))

  # Unavailable parameters warn and are dropped
  expect_warning(partial <- coda::as.mcmc(fit_smooth, pars = c("dhs_phi", "nope")),
                 "not available in model output")
  expect_equal(colnames(partial), "dhs_phi")

  # No available parameters is an error
  expect_error(coda::as.mcmc(fit_smooth, pars = "nope"),
               "None of the requested parameters")
})


test_that("as.mcmc.dsp rejects objects that are not dsp fits", {
  expect_error(as.mcmc.dsp(list(a = 1)), "must be an object of class 'dsp'")
})


test_that("coda diagnostics accept the coerced object", {

  draws <- coda::as.mcmc(fit_smooth, pars = c("dhs_phi", "dhs_mean"))

  expect_length(coda::effectiveSize(draws), 2)
  expect_s3_class(window(draws, start = stats::start(draws) + 100), "mcmc")
})
