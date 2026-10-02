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
