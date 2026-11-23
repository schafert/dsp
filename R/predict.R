#----------------------------------------------------------------------------
#' Posterior predictive summaries
#'
#' @param object object of class dsp from [dsp_fit()]
#' @param cp_thres,cp_prop Deprecated. Changepoint
#'   detection has moved to [changepoints()]. If either argument is supplied
#'   for a changepoint model, `predict()` warns and returns the detected
#'   changepoints as in earlier versions; see Value.
#' @param ... `probs` may be supplied to choose the quantiles returned by
#'   [fitted.dsp()].
#' @details
#' For every model, this method returns `fitted(object, pars = "ypred", ...)`:
#' pointwise summaries of the posterior predictive distribution of a new
#' observation at each fitted time point. Its quantiles give prediction
#' intervals, which include the observation noise. For summaries of the
#' conditional mean \eqn{\mu_t} instead, use `fitted(object, pars = "mu")`.
#'
#' Missing values in `y` are drawn from the posterior predictive distribution
#' during fitting, so `predict()` imputes them. In particular, padding `y`
#' with `NA` at the end before fitting gives forecasts at those time points.
#' For regression, the covariates must be supplied in `X` for every time
#' point, including those with missing `y`; prediction at new covariate
#' values without refitting is not supported.
#'
#' Changepoints of a changepoint model are returned by [changepoints()].
#'
#' @returns
#' A matrix with one row per time point of the fitted series, containing the
#' posterior predictive mean and quantiles. The default quantiles are 2.5%,
#' 25%, 50%, 75%, and 97.5%.
#'
#' If the deprecated `cp_thres` or `cp_prop` is supplied for a changepoint
#' model, the return value of earlier versions: if `cp_prop = FALSE`, an
#' integer vector of the indices of detected changepoints; if
#' `cp_prop = TRUE`, a list with these indices in `cp_t` and, in `cp_prop`,
#' the posterior probabilities of threshold exceedance at all `T` time
#' points, with `NA` at the first `D`.
#'
#' @examples
#' set.seed(200)
#' signal = c(rep(0, 50), rep(10, 50))
#' y = signal + rnorm(100)
#'
#' # Withhold the last five observations to forecast them
#' y_f <- y
#' y_f[96:100] <- NA
#' spec_smooth <- dsp_spec(family = "gaussian", model = "smoothing", D = 1)
#' fit_smooth <- dsp_fit(y_f, spec_smooth, nsave = 200, nburn = 200)
#' tail(predict(fit_smooth, probs = c(0.025, 0.975)), 5)
#'
#' # The conditional mean, without observation noise
#' tail(fitted(fit_smooth, pars = "mu", probs = c(0.025, 0.975)), 5)
#'
#' @seealso [fitted.dsp()], [changepoints()]
#' @method predict dsp
#' @export
predict.dsp <- function(object, cp_thres = deprecated(), cp_prop = deprecated(), ...){

  check_dsp(object)

  old_args <- lifecycle::is_present(cp_thres) || lifecycle::is_present(cp_prop)
  if(old_args && object$model_spec$model == "changepoint"){
    arg <- if(lifecycle::is_present(cp_thres)) "cp_thres" else "cp_prop"
    lifecycle::deprecate_warn(
      "1.8.0", sprintf("predict.dsp(%s)", arg), "changepoints()",
      details = "`predict()` now returns posterior summaries of the mean for every model.")
    if(!lifecycle::is_present(cp_thres)) cp_thres <- 0.5
    if(!lifecycle::is_present(cp_prop)) cp_prop <- FALSE
    probability <- cp_probability(object, cp_thres)
    cp_t <- which(!is.na(probability) & probability >= cp_thres)
    if(isTRUE(cp_prop)) return(list(cp_t = cp_t, cp_prop = probability))
    return(cp_t)
  }
  if(old_args)
    stop("cp_thres and cp_prop apply only to changepoint models; see changepoints().",
         call. = FALSE)

  dots <- list(...)
  if (length(dots) &&
      (is.null(names(dots)) || any(!names(dots) %in% "probs")))
    stop("Only 'probs' is supported through ...; to forecast, set the future values of y to NA before fitting.",
         call. = FALSE)
  if (!"ypred" %in% par_names(object))
    stop("predict() requires the posterior predictive draws 'ypred', which were not saved for this fit.",
         call. = FALSE)
  fitted(object, pars = "ypred", ...)
}
