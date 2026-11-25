#----------------------------------------------------------------------------
#' Pointwise credible intervals
#'
#' Computes equal-tailed pointwise credible intervals for a stored parameter.
#'
#' @param object an object of class `dsp`; see [dsp-class]
#' @param pars character string naming one stored parameter; see
#'   [par_names()]. Defaults to \code{"mu"}.
#' @param level the posterior probability covered by each interval. Defaults
#'   to 0.95.
#' @param times optional vector of length \code{T} labelling the time points
#'   of the series, such as dates. Defaults to \code{1:T}.
#' @param ... currently not used
#'
#' @details
#' Each interval runs from the \code{(1 - level)/2} to the
#' \code{(1 + level)/2} posterior quantile of the parameter at one time point,
#' so it covers that value with posterior probability \code{level}. The
#' intervals are computed at each time point separately. The simultaneous
#' bands of [credible_bands()] instead contain the whole trajectory with
#' posterior probability \code{level}, and are therefore wider.
#'
#' @returns
#' A data frame with columns \code{mean}, \code{lower} and \code{upper}, the
#' posterior mean and the interval limits. For a time-varying parameter there
#' is one row per time point, and a column \code{time}, taken from
#' \code{times}, comes first. For a parameter with one value per predictor and
#' time point, such as the regression coefficients \code{beta}, a column
#' \code{predictor} follows \code{time}, naming the predictor as in
#' [coef.dsp()]. The increments \code{omega} have
#' \code{T - D} entries, which are labelled with the last \code{T - D} time
#' points. The basis coefficients of a B-spline model, \code{beta} and
#' \code{evol_sigma_t2}, are labelled by a column \code{coefficient}, and a
#' scalar parameter gives a single row.
#'
#' @seealso [credible_bands()] for simultaneous bands, [fitted.dsp()] for
#'   other posterior quantiles, [predict.dsp()] for prediction intervals
#'
#' @examples
#' set.seed(200)
#' y <- c(rep(0, 50), rep(10, 50)) + rnorm(100)
#' model_spec <- dsp_spec(family = "gaussian", model = "changepoint", D = 1)
#' fit <- dsp_fit(y, model_spec = model_spec, nsave = 200, nburn = 200)
#'
#' # 95% pointwise intervals for the mean function
#' head(credible_intervals(fit))
#'
#' # 80% intervals for the increments, labelled by date
#' dates <- seq(as.Date("2000-01-01"), by = "month", length.out = 100)
#' head(credible_intervals(fit, pars = "omega", level = 0.8, times = dates))
#'
#' @export
credible_intervals <- function(object, pars = "mu", level = 0.95, times = NULL, ...){

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
  if(length(d) > 3)
    stop("Unsupported parameter shape: must be 2D or 3D.", call. = FALSE)

  probs <- c((1 - level) / 2, (1 + level) / 2)
  intervals_one <- function(b){
    q <- apply(b, 2, stats::quantile, probs = probs, names = FALSE)
    data.frame(mean = colMeans(b), lower = q[1, ], upper = q[2, ])
  }

  if(is.null(d) || length(d) == 1)
    return(intervals_one(matrix(samps, ncol = 1)))

  index <- param_index(object, pars, d[2], times)

  if(length(d) == 2){
    out <- data.frame(index$value, intervals_one(samps))
  } else {
    predictor <- predictor_names(object, d[3])
    out <- do.call(rbind, lapply(seq_len(d[3]), function(j)
      data.frame(index$value, predictor = predictor[j],
                 intervals_one(matrix(samps[, , j], nrow = d[1])))))
  }
  names(out)[1] <- index$name

  rownames(out) <- NULL
  out

}
