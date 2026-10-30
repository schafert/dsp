#----------------------------------------------------------------------------
#' Simultaneous credible bands and band scores
#'
#' Computes simultaneous credible bands for a time-varying parameter, and the
#' simultaneous band scores (SimBaS) that test, at each time point, whether
#' the parameter is zero.
#'
#' @param object an object of class `dsp`; see [dsp-class]
#' @param pars character string naming one time-varying parameter; see
#'   [par_names()]. Defaults to \code{"mu"}.
#' @param level the credibility level of the bands. Defaults to 0.95.
#' @param ... currently not used
#'
#' @details
#' The bands of Crainiceanu et al. (2007) have the form posterior mean plus or
#' minus a multiple of the pointwise posterior standard deviation. The
#' multiple is the \code{level} quantile, over the posterior draws, of the
#' largest standardized absolute deviation from the posterior mean across all
#' time points. The bands therefore contain the whole trajectory in a
#' proportion \code{level} of the posterior draws, whereas the pointwise
#' intervals from [fitted.dsp()] are computed at each time point separately.
#'
#' The SimBaS of Meyer et al. (2015) at time \eqn{t} is the smallest
#' \eqn{\alpha} for which the \eqn{1 - \alpha} simultaneous band excludes zero
#' at \eqn{t}. A small value is evidence that the parameter is nonzero at
#' \eqn{t}, adjusted for testing all time points at once. The minimum over
#' time, \code{min(bands[, "simBaS"], na.rm = TRUE)}, is the global Bayesian
#' p-value for the hypothesis that the parameter is zero at every time point.
#'
#' The scores are informative only for a parameter for which zero is a value
#' of interest. Examples are the increments \code{omega}, where a small score
#' marks a time at which the mean function changes (for \code{D = 1}) or its
#' slope changes (for \code{D = 2}), and the regression coefficients
#' \code{beta}, where a small score marks a time at which the predictor has an
#' effect. For the mean function \code{mu} of a series far from zero, every
#' score is close to zero and carries no information. The changepoint model
#' stores \code{omega} by default. For \code{omega}, the test is distinct
#' from the changepoint rule of [predict.dsp()], which compares the
#' increments with the threshold of the changepoint model.
#'
#' Time points with no posterior variation, such as fixed initial values, are
#' excluded when the bands are calibrated. At those points both bands equal
#' the constant value and the score is \code{NA}.
#'
#' @returns
#' For a time-varying parameter, a numeric matrix with one row per time point
#' (or index) and columns \code{mean}, \code{lower}, \code{upper} and
#' \code{simBaS}. A parameter with a further dimension, such as the
#' regression coefficients \code{beta}, gives an array whose third dimension
#' holds these columns, with one simultaneous band per predictor.
#'
#' @references
#' Crainiceanu, C. M., Ruppert, D., Carroll, R. J., Joshi, A. and Goodner, B.
#' (2007). Spatially adaptive Bayesian penalized splines with heteroscedastic
#' errors. \emph{Journal of Computational and Graphical Statistics}, 16(2),
#' 265--288. \doi{10.1198/106186007X208768}
#'
#' Meyer, M. J., Coull, B. A., Versace, F., Cinciripini, P. and Morris, J. S.
#' (2015). Bayesian function-on-function regression for multilevel functional
#' data. \emph{Biometrics}, 71(3), 563--574. \doi{10.1111/biom.12299}
#'
#' @seealso [fitted.dsp()] for pointwise intervals, [plot.dsp()] with
#'   \code{include_joint_bands = TRUE} to display the bands
#'
#' @examples
#' set.seed(200)
#' y <- c(rep(0, 50), rep(10, 50)) + rnorm(100)
#' model_spec <- dsp_spec(family = "gaussian", model = "changepoint", D = 1)
#' fit <- dsp_fit(y, model_spec = model_spec, nsave = 200, nburn = 200)
#'
#' bands <- credible_bands(fit, pars = "omega")
#' head(bands)
#'
#' # global Bayesian p-value for no change anywhere
#' min(bands[, "simBaS"], na.rm = TRUE)
#'
#' # increments that are nonzero at the joint 5% level
#' which(bands[, "simBaS"] < 0.05)
#'
#' @export
credible_bands <- function(object, pars = "mu", level = 0.95, ...){

  check_dsp(object)

  if(!is.character(pars) || length(pars) != 1)
    stop("'pars' must be a single parameter name.", call. = FALSE)

  if(!is.numeric(level) || length(level) != 1 || is.na(level) ||
     level <= 0 || level >= 1)
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)

  draws <- get_draws(object)

  if(!pars %in% names(draws))
    stop("'pars' must name one of the parameters in the fitted model: ",
         paste(names(draws), collapse = ", "), ".", call. = FALSE)

  samps <- draws[[pars]]
  d <- dim(samps)

  if(is.null(d) || length(d) == 1)
    stop("'", pars, "' is a scalar parameter; credible bands require a ",
         "time-varying parameter.", call. = FALSE)

  bands_one <- function(b){
    cb <- credBands(b, alpha = 1 - level)
    cbind(mean = colMeans(b), lower = cb[, 1], upper = cb[, 2],
          simBaS = simBaS(b))
  }

  if(length(d) == 2)
    return(bands_one(samps))

  if(length(d) == 3){
    out <- array(NA_real_, dim = c(d[2], d[3], 4))
    for(j in seq_len(d[3]))
      out[, j, ] <- bands_one(matrix(samps[, , j], nrow = d[1]))
    dimnames(out) <- list(NULL, NULL, c("mean", "lower", "upper", "simBaS"))
    return(out)
  }

  stop("Unsupported parameter shape: must be 2D or 3D.", call. = FALSE)

}
