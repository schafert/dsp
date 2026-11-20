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
#' @param times optional vector of length \code{T} labelling the time points
#'   of the series, such as dates. Defaults to \code{1:T}.
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
#' An object of class \code{dsp_bands}, a data frame with one row per time
#' point and columns \code{time}, \code{mean}, \code{lower}, \code{upper} and
#' \code{simBaS}. For a parameter with one band per predictor, such as the
#' regression coefficients \code{beta}, a column \code{predictor} follows
#' \code{time}. The \code{time} column is taken from \code{times}. The
#' increments \code{omega} have \code{T - D} entries, which are labelled with
#' the last \code{T - D} time points, so that each row refers to the time at
#' which the increment ends. The print method reports the global Bayesian
#' p-value and the periods in which the band excludes zero, by predictor.
#' Subsetting with \code{[} returns a plain data frame.
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
#' # global Bayesian p-value and the periods in which the band excludes zero
#' bands <- credible_bands(fit, pars = "omega")
#' bands
#'
#' # the pointwise values
#' head(as.data.frame(bands))
#'
#' # time points at which the increment is nonzero at the joint 5% level
#' subset(bands, simBaS < 0.05)$time
#'
#' @export
credible_bands <- function(object, pars = "mu", level = 0.95, times = NULL, ...){

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
  if(length(d) > 3)
    stop("Unsupported parameter shape: must be 2D or 3D.", call. = FALSE)

  nT <- nobs(object)
  if(is.null(times)) times <- seq_len(nT)
  if(length(times) != nT)
    stop("'times' must have one entry per time point of the series (", nT, ").",
         call. = FALSE)

  # A parameter with fewer than T entries, such as the increments omega with
  # T - D, is aligned with the last time points
  n <- d[2]
  if(n > nT)
    stop("'", pars, "' has more entries than the series has time points.",
         call. = FALSE)
  time <- times[(nT - n + 1):nT]

  bands_one <- function(b){
    cb <- credBands(b, alpha = 1 - level)
    data.frame(mean = colMeans(b), lower = cb[, 1], upper = cb[, 2],
               simBaS = simBaS(b))
  }

  if(length(d) == 2){
    out <- data.frame(time = time, bands_one(samps))
  } else {
    X <- get_spec(object)$arguments$X
    predictor <- if(!is.null(colnames(X))) colnames(X) else seq_len(d[3])
    out <- do.call(rbind, lapply(seq_len(d[3]), function(j)
      data.frame(time = time, predictor = predictor[j],
                 bands_one(matrix(samps[, , j], nrow = d[1])))))
  }

  rownames(out) <- NULL
  structure(out, class = c("dsp_bands", "data.frame"),
            pars = pars, level = level, ndraws = d[1])

}

#' @rdname credible_bands
#' @param x an object of class \code{dsp_bands}
#' @param digits number of significant digits for the global Bayesian p-value
#' @param max_periods the largest number of periods to list per band
#' @method print dsp_bands
#' @export
print.dsp_bands <- function(x, digits = 2, max_periods = 10, ...){

  level <- attr(x, "level")
  df <- as.data.frame(x)
  groups <- if(is.null(df$predictor)) list(df) else split(df, factor(df$predictor, levels = unique(df$predictor)))

  cat("Simultaneous ", format(100*level), "% credible bands for ", attr(x, "pars"),
      " (", nrow(groups[[1]]), " time points)\n", sep = "")

  for(g in seq_along(groups)){
    b <- groups[[g]]
    if(!is.null(df$predictor)) cat("\nPredictor ", names(groups)[g], ":\n", sep = "")
    p <- suppressWarnings(min(b$simBaS, na.rm = TRUE))
    # The scores are resolved to 1/ndraws, so a zero score is reported as a bound
    pval <- if(!is.finite(p)) "NA" else if(p == 0)
      paste("<", format(signif(1/attr(x, "ndraws"), digits))) else format(signif(p, digits))
    cat("Global Bayesian p-value (minimum SimBaS):", pval, "\n")

    excl <- !is.na(b$simBaS) & (b$lower > 0 | b$upper < 0)
    if(!any(excl)){
      cat("The band includes zero at every time point.\n")
      next
    }
    runs <- rle(excl)
    end <- cumsum(runs$lengths)
    start <- end - runs$lengths + 1
    periods <- data.frame(start = b$time[start[runs$values]],
                          end = b$time[end[runs$values]],
                          length = runs$lengths[runs$values])
    cat("The band excludes zero at ", sum(excl), " time point",
        if(sum(excl) > 1) "s", ", in ", nrow(periods), " period",
        if(nrow(periods) > 1) "s", ":\n", sep = "")
    print(utils::head(periods, max_periods), row.names = FALSE)
    if(nrow(periods) > max_periods)
      cat("... and", nrow(periods) - max_periods, "more periods\n")
  }

  cat("\nUse as.data.frame() for the values at each time point.\n")
  invisible(x)
}

#' @export
`[.dsp_bands` <- function(x, ...){
  out <- NextMethod()
  if(is.data.frame(out)){
    class(out) <- "data.frame"
    attr(out, "pars") <- NULL
    attr(out, "level") <- NULL
    attr(out, "ndraws") <- NULL
  }
  out
}

#' @export
as.data.frame.dsp_bands <- function(x, ...){
  class(x) <- "data.frame"
  attr(x, "pars") <- NULL
  attr(x, "level") <- NULL
  attr(x, "ndraws") <- NULL
  x
}
