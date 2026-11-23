#' Extract posterior coefficient estimates
#'
#' Return posterior means of the regression or B-spline coefficients.
#'
#' @param object a fitted `dsp` object.
#' @param ... currently unused.
#' @returns For regression, a matrix with one row per fitted time point and
#'   one column per predictor, named by the column names of `X`, or `x1`,
#'   `x2`, ... if `X` has none. For B-splines, a named vector of basis
#'   coefficient estimates. Other model types do not have a separate
#'   coefficient parameter; use [fitted.dsp()] for their mean process.
#' @seealso [fitted.dsp()], [confint.dsp()]
#' @importFrom stats coef
#' @method coef dsp
#' @export
coef.dsp <- function(object, ...) {
  check_dsp(object)
  if (!object$model_spec$model %in% c("regression", "bspline"))
    stop("Coefficient extraction is available for regression and B-spline models; use fitted() for the mean process.",
         call. = FALSE)

  beta <- get_draws(object)[["beta"]]
  if (is.null(beta))
    stop("Posterior draws of 'beta' were not saved.", call. = FALSE)

  d <- dim(beta)
  if (object$model_spec$model == "regression") {
    if (length(d) != 3L)
      stop("Regression coefficient draws must have draw, time, and predictor dimensions.",
           call. = FALSE)
    out <- matrix(colMeans(matrix(beta, nrow = d[1L])),
                  nrow = d[2L], ncol = d[3L])
    X <- object$model_spec$arguments$X
    predictor_names <- colnames(X)
    if (is.null(predictor_names))
      predictor_names <- paste0("x", seq_len(d[3L]))
    dimnames(out) <- list(rownames(X), predictor_names)
    return(out)
  }

  if (length(d) != 2L)
    stop("B-spline coefficient draws must be a draw-by-coefficient matrix.",
         call. = FALSE)
  out <- colMeans(beta)
  names(out) <- get_flatnames("beta", d[-1L])
  out
}

#' Extract response residuals
#'
#' Subtract the posterior mean of the stored `mu` process from the original
#' observations. For changepoint models, these residuals include any additive
#' anomaly component rather than subtracting it separately.
#'
#' @param object a fitted `dsp` object.
#' @param y original observations. Defaults to the observations retained by
#'   [dsp_fit()]. Supply this argument for older fits that did not store them.
#' @param ... currently unused.
#' @returns A vector of response residuals, with `NA` at originally missing
#'   observations. No standardization by the observation variance is applied.
#' @seealso [fitted.dsp()]
#' @importFrom stats residuals
#' @method residuals dsp
#' @export
residuals.dsp <- function(object, y = object$y, ...) {
  check_dsp(object)
  if (is.null(y))
    stop("Original observations were not stored in this fit; supply them with y or refit with the updated dsp_fit().",
         call. = FALSE)
  mu <- get_draws(object)[["mu"]]
  if (is.null(mu) || length(dim(mu)) != 2L)
    stop("Posterior draws of 'mu' must be saved as a draw-by-time matrix.",
         call. = FALSE)
  if (!is.numeric(y) || length(y) != ncol(mu))
    stop("'y' must contain one numeric observation per fitted time point.",
         call. = FALSE)
  y - colMeans(mu)
}

#' Extract posterior log-likelihood values
#'
#' Return the posterior mean log-likelihood or the stored log-likelihood
#' value at each saved MCMC iteration.
#'
#' @param object a fitted `dsp` object.
#' @param summary either `"mean"` (the default) or `"draws"`.
#' @param ... currently unused.
#' @returns For `summary = "mean"`, an object of class `logLik` with
#'   attributes `nobs` and `df = NA`. For `summary = "draws"`, a numeric
#'   vector with one value per saved iteration.
#' @details The mean is the average of the stored log-likelihood values,
#'   not the log-likelihood evaluated at posterior mean parameters or the
#'   logarithm of the average likelihood. It is not a maximized
#'   log-likelihood, and a shrinkage model has no fixed parameter count, so
#'   the degrees of freedom are `NA` and `AIC()` and `BIC()` return `NA`.
#'   Use [DIC()] for the model-comparison summaries provided by the package.
#' @seealso [as.mcmc.dsp()], [DIC()], [deviance.dsp()]
#' @importFrom stats logLik
#' @method logLik dsp
#' @export
logLik.dsp <- function(object, summary = c("mean", "draws"), ...) {
  check_dsp(object)
  summary <- match.arg(summary)
  draws <- get_draws(object)[["loglike"]]
  if (is.null(draws))
    stop("Posterior log-likelihood draws were not saved.", call. = FALSE)
  draws <- as.numeric(draws)
  if (summary == "draws") return(draws)
  structure(mean(draws), df = NA_real_, nobs = nobs(object), class = "logLik")
}

#' Extract the posterior mean deviance
#'
#' Return the posterior mean of the deviance, \eqn{\bar{D} = -2} times the
#' average of the stored log-likelihood values.
#'
#' @param object a fitted `dsp` object.
#' @param ... currently unused.
#' @returns A numeric scalar.
#' @details \eqn{\bar{D}} is the goodness-of-fit term of the deviance
#'   information criterion: with the \eqn{p_D} estimator of the effective
#'   number of parameters, \eqn{DIC = \bar{D} + p_D}. [DIC()] returns both
#'   versions of the criterion.
#' @seealso [DIC()], [logLik.dsp()]
#' @importFrom stats deviance
#' @method deviance dsp
#' @export
deviance.dsp <- function(object, ...) {
  check_dsp(object)
  -2 * as.numeric(logLik(object))
}

#' Extract pointwise posterior credible intervals
#'
#' The `confint()` method for `dsp` returns equal-tailed Bayesian credible
#' intervals, not frequentist confidence intervals or simultaneous bands.
#'
#' @param object a fitted `dsp` object.
#' @param parm one stored parameter name from [par_names()]. Defaults to
#'   `"mu"`, as in [fitted.dsp()].
#' @param level posterior probability covered by each interval. Defaults to
#'   `0.95`.
#' @param ... `pars` may be used in place of `parm`, for consistency with
#'   [fitted.dsp()] and [as.mcmc.dsp()]. Other arguments are not accepted.
#' @returns A matrix with one row per selected parameter element and two
#'   columns for the lower and upper posterior quantiles.
#' @seealso [fitted.dsp()], [credible_bands()], [coef.dsp()]
#' @importFrom stats confint
#' @method confint dsp
#' @export
confint.dsp <- function(object, parm = "mu", level = 0.95, ...) {
  check_dsp(object)
  dots <- list(...)
  if ("pars" %in% names(dots)) {
    if (!missing(parm))
      stop("Supply the parameter through 'parm' or 'pars', not both.", call. = FALSE)
    parm <- dots$pars
    dots$pars <- NULL
  }
  if (length(dots))
    stop("Unused argument(s) to confint(): ",
         paste(if (is.null(names(dots))) "unnamed" else names(dots), collapse = ", "),
         call. = FALSE)
  pars <- parm
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1)
    stop("'level' must be a number strictly between zero and one.",
         call. = FALSE)
  if (!is.character(pars) || length(pars) != 1L || is.na(pars))
    stop("'parm' must be a single stored parameter name.", call. = FALSE)
  draws <- as.matrix(as.mcmc(object, pars = pars))

  probs <- c((1 - level) / 2, (1 + level) / 2)
  out <- t(vapply(seq_len(ncol(draws)), function(j) {
    stats::quantile(draws[, j], probs = probs, names = FALSE)
  }, numeric(2L)))
  dimnames(out) <- list(colnames(draws), names(stats::quantile(0, probs)))
  out
}
