#----------------------------------------------------------------------------
#' Predict changepoints from the output of ABCO
#'
#' @param object object of class dsp from [dsp_fit()]
#' @param cp_thres (default 0.5) cutoff proportion for percentage of posterior
#'                  samples exceeding the threshold needed to label a changepoint
#' @param cp_prop (default FALSE) logical flag determining if the posterior proportions of threshold exceedance is to be returned.
#' @param ... currently unused
#' @details
#' The changepoint model uses a thresholding mechanism with a latent indicator variable.
#' This function calculates the proportion of samples where the increment exceeds the threshold.
#'
#'
#' @returns
#' If cp_prop = FALSE, a numeric vector of indices that correspond to indices of the observed data.
#' If cp_prop = TRUE, a list containing:
#'
#'     - 'cp_t':  a numeric vector of indices that correspond to indices of the observed data.
#'     - 'cp_prop': a numeric vector of length (T - D) with the pointwise proportion of samples where the increment exceeds the threshold.
#'
#' If no proportions exceed cp_thres, then the vector will be a length 0 integer vector.
#'
#' @examples
#'
#' signal = c(rep(0, 50), rep(10, 50))
#' noise = rep(1, 100)
#' noise_var = rep(1, 100)
#' for (k in 2:100){
#'   noise_var[k] = exp(0.9*log(noise_var[k-1]) + rnorm(1, 0, 0.5))
#'   noise[k] = rnorm(1, 0, sqrt(noise_var[k])) }
#'
#' y = signal + noise
#' model_spec = dsp_spec(family = "gaussian", model = "changepoint")
#' mcmc_output = dsp_fit(y, model_spec = model_spec)
#' predict(mcmc_output)
#'
#'
#' @export

predict.dsp <- function(object, cp_thres = 0.5, cp_prop = FALSE, ...){

  if(object$model_spec$model != 'changepoint'){
    stop("Changepoint model was not fit.")
  }

  if(!all(c("omega", "r") %in% names(object$mcmc_output))){
    stop("Can't calculate changepoint evidence without samples of omega and r.")
  }

  if (is.na(cp_thres) || !is.numeric(cp_thres) || cp_thres < 0 || cp_thres > 1) stop("cp_thres must be a numeric value between 0 and 1.")

  cp_list = rep(0, length(object$mcmc_output$omega[1,]))
  D = object$model_spec$arguments$D

  cp_mat <- apply(object$mcmc_output$omega^2, MARGIN = 2, \(samp){samp > exp(object$mcmc_output$r)})

  cp_list <- colMeans(cp_mat)

  cp_t <- which(cp_list >= cp_thres) + D

  if(cp_prop){
    return(list(cp_t = cp_t, cp_prop = cp_list))
  }

  return(cp_t)

}

