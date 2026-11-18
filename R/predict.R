#----------------------------------------------------------------------------
#' Extract in-sample predictions or detect changepoints
#'
#' @param object object of class dsp from [dsp_fit()]
#' @param cp_thres (default 0.5) cutoff proportion for percentage of posterior
#'                  samples exceeding the threshold needed to label a changepoint.
#'                  Used only for changepoint models.
#' @param cp_prop (default FALSE) logical flag determining if the posterior proportions of threshold exceedance is to be returned.
#'   Used only for changepoint models.
#' @param ... for regression, smoothing, and B-spline models, `probs` may be
#'   supplied to choose the posterior quantiles returned by [fitted.dsp()].
#'   Currently unused for changepoint models.
#' @details
#' For regression, smoothing (Gaussian or negative binomial), and B-spline
#' models, this method returns `fitted(object, pars = "mu", ...)`: posterior
#' summaries of the conditional mean at the fitted observation times. For
#' regression, these use the original covariates and the coefficient draws
#' at their corresponding fitted time points. New covariates, new observation
#' times, and forecasts beyond the fitted series are not supported.
#' These summaries describe the conditional mean, not new observation draws;
#' use `fitted(object, pars = "ypred")` for posterior predictive summaries.
#'
#' The changepoint model uses a thresholding mechanism with a latent indicator variable.
#' For this model, the method instead calculates the proportion of samples
#' where the increment exceeds the threshold and returns detected changepoints.
#'
#'
#' @returns
#' For regression, smoothing, and B-spline models, a matrix with one row per
#' fitted time point, containing the posterior mean and quantiles of `mu`.
#' The default quantiles are 2.5%, 25%, 50%, 75%, and 97.5%.
#'
#' For changepoint models, if `cp_prop = FALSE`, an integer vector containing the indices of detected changepoints in the observed data.
#' If no changepoints are detected, an empty integer vector is returned.
#'
#' If cp_prop = TRUE, a list containing:
#'     - 'cp_t': a data frame with 'index' and 'probability' columns for detected changepoints.
#'     - 'cp_prop': a numeric vector of length `T` containing the pointwise posterior proportions of threshold exceedance,
#'       aligned with the observed series.
#'       The first `D` entries are `NA`, because no increment is defined there.
#'
#' @examples
#' set.seed(200)
#' signal = c(rep(0, 50), rep(10, 50))
#' noise = rep(1, 100)
#' noise_var = rep(1, 100)
#' for (k in 2:100){
#'   noise_var[k] = exp(0.9*log(noise_var[k-1]) + rnorm(1, 0, 0.5))
#'   noise[k] = rnorm(1, 0, sqrt(noise_var[k])) }
#'
#' y = signal + noise
#' model_spec = dsp_spec(family = "gaussian", model = "changepoint",
#'                       D = 1, useAnom = TRUE)
#' mcmc_output = dsp_fit(y, model_spec = model_spec, nsave = 200, nburn = 200)
#' predict(mcmc_output)
#'
#' predict(mcmc_output, cp_prop = TRUE)
#'
#' # In-sample conditional-mean summaries for a smoothing model
#' spec_smooth <- dsp_spec(family = "gaussian", model = "smoothing", D = 1)
#' fit_smooth <- dsp_fit(y, spec_smooth, nsave = 200, nburn = 200)
#' head(predict(fit_smooth))
#' head(predict(fit_smooth, probs = c(0.025, 0.975)))
#'
#' @seealso [fitted.dsp()]
#' @method predict dsp
#' @export

predict.dsp <- function(object, cp_thres = 0.5, cp_prop = FALSE, ...){

  check_dsp(object)

  if(object$model_spec$model != 'changepoint'){
    dots <- list(...)
    if (length(dots) &&
        (is.null(names(dots)) || any(!names(dots) %in% "probs")))
      stop("Only 'probs' is supported through ... for in-sample predictions; new data and forecasting are not supported.",
           call. = FALSE)
    return(fitted(object, pars = "mu", ...))
  }

  if(!all(c("omega", "gamma") %in% names(object$mcmc_output))){
    stop("Can't calculate changepoint evidence without samples of omega and gamma.")
  }

  if (is.na(cp_thres) || !is.numeric(cp_thres) || cp_thres < 0 || cp_thres > 1) stop("cp_thres must be a numeric value between 0 and 1.")

  D = object$model_spec$arguments$D

  cp_mat <- apply(object$mcmc_output$omega^2, MARGIN = 2, \(samp){samp > exp(object$mcmc_output$gamma)})

  cp_list <- colMeans(cp_mat)

  cp_t <- as.integer(which(cp_list >= cp_thres) + D)

  # if(cp_prop){
  #   # omega holds the T - D increments, so pad the front to index by time
  #   return(list(cp_t = cp_t, cp_prop = c(rep(NA_real_, D), cp_list)))
  # }
  if (cp_prop) {
  probability <- c(rep(NA_real_, D), cp_list)

  return(list(
    cp_t = data.frame(
      index = cp_t,
      probability = probability[cp_t]
    ),
    cp_prop = probability
  ))
}

  return(cp_t)

}
