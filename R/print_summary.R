print.dsp <- function(object, ...){



}


#----------------------------------------------------------------------------
#' Predict changepoints from the output of ABCO
#'
#' @param object object of class 'dsp'
#' @param cp_thres (default 0.5) cutoff proportion for percentage of posterior
#'                  samples exceeding the threshold needed to label a changepoint
#' @param cp_prop (default FALSE) logical flag determining if the posterior proportions of threshold exceedance is to be returned.
#'
#' @details
#' The changepoint model uses a thresholding mechanism with a latent indicator variable.
#' This function calculates the proportion of samples where the increment exceeds the threshold.
#'
#'
#' @return
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

  if (is.na(cp_thres) || !is.numeric(cp_thres) || cp_thres < 0 || cp_thres > 1) stop("cp_thres must be a numeric value between 0 and 1.") ,

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

summary.dsp <- function(object, pars, probs = c(0.025, 0.25, 0.50, 0.75, 0.975), ...){

  pars_in <- pars %in% names(object$mcmc_output)

  if(!all(pars_in)){

    warning(paste("The following parameters were requested, but not available in model output:", paste(pars[!pars_in], collapse = ", ")))

    pars <- pars[pars_in]
  }



}

