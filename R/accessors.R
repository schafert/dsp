#----------------------------------------------------------------------------
#' Accessor functions for fitted dsp models
#'
#' Extract information from a fitted \code{dsp} object without referring to
#' its internal structure.
#'
#' @param object an object of class `dsp`; see [dsp-class]
#' @param pars character string naming one parameter; see [par_names()] for
#'   the parameters available in \code{object}. Defaults to \code{"mu"}, the
#'   mean function, which every model returns.
#' @param probs numeric vector of probabilities for the posterior quantiles
#' @param ... currently not used
#'
#' @details
#' \code{par_names()} returns the names of the parameters with stored posterior
#' draws. These are the values accepted by the \code{pars} argument of
#' \code{fitted()}, [summary.dsp()] and [as.mcmc.dsp()], and by the
#' \code{type} argument of [plot.dsp()].
#'
#' \code{fitted()} returns pointwise posterior summaries of one parameter.
#' For the posterior draws themselves, use [as.mcmc.dsp()].
#'
#' \code{DIC()} returns the deviance information criterion under the two
#' standard estimators of the effective number of parameters,
#' \eqn{p_D} (Spiegelhalter et al. 2002) and \eqn{p_V} (Gelman et al. 2014).
#' Lower values are preferred. We recommend \eqn{p_D} by default.
#'
#' \code{get_spec()} returns the [dsp_spec()] object that the model was
#' fitted with.
#'
#' \code{nobs()} returns the number of time points \code{T} in the fitted
#' series, including any missing observations, which the model imputes.
#'
#' @returns
#' \code{par_names()}: a character vector.
#'
#' \code{fitted()}: for a time-varying parameter, a numeric matrix with one row
#' per time point (or index) and columns \code{"mean"} followed by one column
#' per element of \code{probs}, named as by [stats::quantile()]. A parameter
#' with a further dimension, such as the regression coefficients \code{beta},
#' gives an array whose third dimension holds these summaries. A scalar
#' parameter gives a named vector.
#'
#' \code{DIC()}: a data frame with one row per estimator and columns
#' \code{estimator}, \code{DIC} and \code{p_eff}, the effective number of
#' parameters.
#'
#' \code{get_spec()}: an object of class \code{dsp_spec}.
#'
#' \code{nobs()}: an integer.
#'
#' @references
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P. and van der Linde, A.
#' (2002). Bayesian measures of model complexity and fit. \emph{Journal of the
#' Royal Statistical Society: Series B}, 64(4), 583--639.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A. and
#' Rubin, D. B. (2014). \emph{Bayesian Data Analysis}, 3rd edition. CRC Press.
#'
#' @seealso [summary.dsp()], [as.mcmc.dsp()], [plot.dsp()], [predict.dsp()]
#'
#' @examples
#' set.seed(200)
#' y <- c(rep(0, 50), rep(10, 50)) + rnorm(100)
#' model_spec <- dsp_spec(family = "gaussian", model = "smoothing", D = 1)
#' fit <- dsp_fit(y, model_spec = model_spec, nsave = 200, nburn = 200)
#'
#' par_names(fit)
#' head(fitted(fit))
#' DIC(fit)
#' get_spec(fit)
#' nobs(fit)
#'
#' @name dsp-accessors
NULL

#' @rdname dsp-accessors
#' @export
par_names <- function(object){

  check_dsp(object)

  names(get_draws(object))

}

#' @rdname dsp-accessors
#' @importFrom stats fitted
#' @method fitted dsp
#' @export
fitted.dsp <- function(object, pars = "mu",
                       probs = c(0.025, 0.25, 0.50, 0.75, 0.975), ...){

  check_dsp(object)

  if(!is.character(pars) || length(pars) != 1)
    stop("'pars' must be a single parameter name.", call. = FALSE)

  draws <- get_draws(object)

  if(!pars %in% names(draws))
    stop("'pars' must name one of the parameters in the fitted model: ",
         paste(names(draws), collapse = ", "), ".", call. = FALSE)

  summarize_draws(draws[[pars]], probs = probs)

}

#' @rdname dsp-accessors
#' @export
DIC <- function(object, ...){
  UseMethod("DIC")
}

#' @rdname dsp-accessors
#' @method DIC dsp
#' @export
DIC.dsp <- function(object, ...){

  check_dsp(object)

  if(is.null(object$DIC$DIC))
    stop("The DIC was not computed for this fit; refit with computeDIC = TRUE.",
         call. = FALSE)

  # abco(), btf() and the other fitters store the p_D version first and the
  # p_V version second.
  data.frame(estimator = c("p_D", "p_V"),
             DIC = object$DIC$DIC,
             p_eff = object$DIC$p_d)

}

#' @rdname dsp-accessors
#' @export
get_spec <- function(object){

  check_dsp(object)

  object$model_spec

}

#' @rdname dsp-accessors
#' @importFrom stats nobs
#' @method nobs dsp
#' @export
nobs.dsp <- function(object, ...){

  check_dsp(object)

  ncol(object$mcmc_output$mu)

}

#----------------------------------------------------------------------------
#' Stop unless an object is a fitted dsp model
#'
#' @param object any \R object
#' @return \code{object}, invisibly.
#' @keywords internal

check_dsp <- function(object){

  if(!is.dsp(object))
    stop("'object' must be an object of class 'dsp'.", call. = FALSE)

  invisible(object)

}

#----------------------------------------------------------------------------
#' Extract the posterior draws from a fitted dsp model
#'
#' @param object an object of class `dsp`
#' @return The elements of \code{object$mcmc_output} that hold posterior draws,
#'   with the draws in the first dimension.
#' @keywords internal

get_draws <- function(object){

  nsave <- unname(object$mcpar["nsave"])

  # DIC and p_d are posterior summaries, and any element whose leading
  # dimension is not the number of saved draws cannot be a monitored quantity.
  draws <- object$mcmc_output[!names(object$mcmc_output) %in% c("DIC", "p_d")]

  is_draw <- vapply(draws, function(samps){
    d <- dim(samps)
    if(is.null(d)) length(samps) == nsave else d[1] == nsave
  }, logical(1))

  draws[is_draw]

}

#----------------------------------------------------------------------------
#' Pointwise posterior summaries of one parameter
#'
#' @param samps posterior draws of one parameter, with the draws in the first
#'   dimension: a vector, a matrix or a three-dimensional array
#' @param probs numeric vector of probabilities for the posterior quantiles
#' @return A named vector for a scalar parameter, a matrix with one row per
#'   index for a vector parameter, or an array whose third dimension holds the
#'   summaries for a matrix parameter.
#' @keywords internal

summarize_draws <- function(samps, probs){

  summary_fun <- function(col) c(mean = mean(col), stats::quantile(col, probs = probs))

  d <- dim(samps)

  # A one-dimensional array, as btf_nb() returns for r, is a scalar parameter
  if(is.null(d) || length(d) == 1)
    return(summary_fun(as.vector(samps)))

  if(length(d) == 2)
    return(t(apply(samps, 2, summary_fun)))

  if(length(d) == 3){
    out <- array(NA_real_, dim = c(d[2], d[3], length(probs) + 1))
    for(j in seq_len(d[3]))
      out[, j, ] <- t(apply(samps[, , j, drop = FALSE][, , 1], 2, summary_fun))
    dimnames(out) <- list(NULL, NULL, names(summary_fun(samps[, 1, 1])))
    return(out)
  }

  stop("Unsupported parameter shape: must be scalar, 2D, or 3D.", call. = FALSE)

}
