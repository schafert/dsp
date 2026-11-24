# A small fit, shared by the tests below.
fit_plot <- suppressMessages(
  dsp_fit(y = {set.seed(1); c(rep(0, 40), rep(6, 40)) + rnorm(80)},
          model_spec = dsp_spec(family = "gaussian", model = "smoothing", D = 1),
          nsave = 100, nburn = 100, verbose = FALSE)
)

test_that("plot.dsp uses its default type when the argument is omitted", {
  # missing() is TRUE exactly when the caller omits the argument, so a guard
  # written as if(missing(type)) makes the default unreachable.
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_silent(plot(fit_plot))
  expect_identical(formals(dsp:::plot.dsp)$type, "mu")
})

test_that("plot.dsp rejects an unusable type with one informative message", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # NULL makes missing() FALSE, so it must be caught by validating the value
  expect_error(plot(fit_plot, type = NULL), "must name one of")
  expect_error(plot(fit_plot, type = "not_a_parameter"), "must name one of")
  expect_error(plot(fit_plot, type = c("mu", "ypred")), "must name one of")
  expect_error(plot(fit_plot, type = 1), "must name one of")

  # the message names the parameters the object actually holds
  expect_error(plot(fit_plot, type = "nope"), "mu")
})

test_that("plot.dsp accepts a supplied type, named or positional", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_silent(plot(fit_plot, type = "mu"))
  expect_silent(plot(fit_plot, "ypred"))
})

test_that("plot.dsp aligns omega and the changepoint markers with the series", {
  y <- {set.seed(200); c(rep(0, 40), rep(8, 40)) + rnorm(80)}
  fit <- suppressWarnings(suppressMessages(
    dsp_fit(y = y, model_spec = dsp_spec(family = "gaussian",
                                         model = "changepoint", D = 2),
            nsave = 100, nburn = 100, verbose = FALSE, seed = 1)
  ))
  dates <- seq(as.Date("2000-01-01"), by = "month", length.out = 80)
  cp <- changepoints(fit, times = dates)
  expect_gt(nrow(cp), 0L)

  # record where the changepoint lines are drawn
  drawn <- list()
  local_mocked_bindings(abline = function(..., v = NULL) {
    if (!is.null(v)) drawn[[length(drawn) + 1L]] <<- v
    invisible(NULL)
  }, .package = "dsp")
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # omega has T - D columns, but full-length times are accepted
  expect_no_error(plot(fit, "omega", times = dates))
  expect_no_error(plot(fit, "mu", times = dates))
  expect_equal(drawn[[1]], cp$time)
  expect_equal(drawn[[2]], cp$time)

  # without times, both types mark the same positions on [0, 1]
  drawn <- list()
  plot(fit, "omega")
  plot(fit, "mu")
  expect_equal(drawn[[1]], drawn[[2]])
})
