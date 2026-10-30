# dsp 1.8.0

* Added accessor functions so that a fitted model can be inspected without
  referring to the structure of the object: `par_names()` lists the parameters
  with stored draws, `fitted()` returns pointwise posterior summaries of one
  parameter, `DIC()` returns the deviance information criterion under both
  estimators of the effective number of parameters, `get_spec()` returns the
  model specification and `nobs()` returns the number of time points.
* `summary()` now returns an object of class `summary.dsp` with a print method.
  It reports the model specification, the MCMC settings, the DIC and one line
  per parameter, and fits on one screen. Effective sample sizes are computed by
  `coda::effectiveSize()`; pointwise ones for time-varying parameters are
  opt in through `ess = TRUE`. The pointwise posterior summaries that
  `summary()` previously returned as a list are now returned by `fitted()`, with
  quantile columns named as by `quantile()` (`"2.5%"` rather than `"q_2.5"`).
* Added `credible_bands()`, which returns the simultaneous credible bands of a
  time-varying parameter and its simultaneous band scores (SimBaS), the
  multiplicity-adjusted evidence at each time point that the parameter is
  nonzero. Previously the bands were available only on plots and SimBaS not at
  all. The result is a data frame of class `dsp_bands` with one row per time
  point, labelled by the optional `times`. Its print method reports the global
  Bayesian p-value and the periods in which the band excludes zero.
* The simultaneous credible bands drawn by `plot(include_joint_bands = TRUE)`
  were undefined when any time point had no posterior variation, such as a
  fixed initial value. Such points are now excluded when the bands are
  calibrated.
* A constant observation variance (`obsSV = "const"`, and every B-spline
  model) is now stored as one draw per iteration rather than copied to every
  time point. `summary()` reports it with the scalar parameters, and `plot()`
  still draws it over time for comparison with time-varying variances.
* `summary()` reports the first `D` entries of `evol_sigma_t2`, the variances
  of the initial states, on their own line.
* The anomaly component of the changepoint model (`useAnom = TRUE`) now has
  the horseshoe+ prior of Wu et al. (2024), as documented. Earlier versions
  used a horseshoe prior, because the horseshoe+ lines were commented out.
* The anomaly `zeta` and its variance `zeta_sigma_t2` are now returned with
  `T` columns aligned with `y`, rather than `T - D`. The anomaly is defined
  for `t > D` only, as in Wu et al. (2024), so both are zero at the first `D`
  time points.
* `plot()` no longer drops the first `D` entries of `zeta_sigma_t2`, which has
  no initial states, or the first entry of `evol_sigma_t2` when `D = 0`.
* `as.mcmc()` is re-exported from `coda`, so `as.mcmc(fit)` works after
  `library(dsp)` alone.
* `predict(cp_prop = TRUE)` now returns the changepoint probabilities aligned
  with the observed series, with `NA` for the first `D` time points, so that
  `cp_prop[cp_t]` gives the probabilities at the detected changepoints.
* `print()` for fitted models reported `nskip` as the thinning interval. It now
  reports `nskip + 1`, the interval `coda::thin()` gives.
* `plot()` now rejects `type = "DIC"` and `type = "p_d"`, which are not draws.
* Fixed a missing "and" in the printed specification of B-spline models.

# dsp 1.7.0

* Added a `coda::as.mcmc()` method for fitted `dsp` objects, so that the
  convergence diagnostics and plotting functions in `coda` can be applied
  directly. Parameters are flattened to one column per index, named so that
  individual time points can be addressed, and the result records the
  sampler's iteration index, including burn-in and thinning.
* Added `is.dsp()` and documented the `dsp` class, including the components of
  a fitted object and their dimensions.
* `plot()` now works without an explicit `type`, defaulting to `"mu"`. The
  default was previously unreachable because the method tested whether the
  argument had been supplied rather than testing its value. Supplying `NULL`,
  an unrecognised name, a vector or a non-character value now gives one message
  naming the quantities the fitted object holds.
* Removed `RemoteType`, `RemotePkgRef` and `RemoteUrl` from DESCRIPTION, which
  were left behind by a local installation.
* Removed two unused internal functions, `getEffSize()` and `ergMean()`.

# dsp 1.6.0

* Fixed random number generation so that `set.seed` controls results; `btf_nb` reset the seed and the state sampler drew from a separate stream.
* Added a `seed` argument to `dsp_fit`.
* Fixed `r_user`, which did not hold the overdispersion parameter fixed as documented.
* Negative binomial models now use the faster Polya-Gamma sampler when the overdispersion is integer valued.
* `dsp_spec` no longer accepts `r_init` and `r_sample` directly; both are set through `r_user`.
* Failed Cholesky factorizations are now detected and redrawn with a warning, rather than returning invalid states.
* Fixed plot method dropping every changepoint after the first and drawing them `D` observations early.
* Fixed pluralization of "degree" in the print method.
* Documented the `slice` and `mh` overdispersion samplers as experimental.

# dsp 1.5.1

* minor changes to colors used in plot for consistency

# dsp 1.5.0

* `dsp_spec` now abstracts overdispersion sampling from the user by simplifying to options to default or Poisson approximation.

# dsp 1.4.1

* Changed typo in print method for family "negbinomial"

# dsp 1.4.0

* Changed argument name t01 to times for consistency
* Changed variable name yhat to ypred for correctness

# dsp 1.3.0

* Fixed inconsistencies in model specification family names.
* Enhanced the plotting method.

# dsp 1.2.0

# dsp 1.1.0

# dsp 1.0.0

* Initial CRAN submission.
