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
