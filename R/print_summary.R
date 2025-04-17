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

#----------------------------------------------------------------------------
#' Summarize DSP MCMC chains
#'
#' @param object object of class dsp from [dsp_fit()]
#' @param pars parameter names specified for summaries; currently defaults to all parameters named in object$mcmc_output
#' @param probs numeric vector of [quantile()]s requested for posterior summary of pars. Defaults to c(0.025, 0.25, 0.50, 0.75, 0.975)
#' @param ... currently unused
#'
#' @returns
#' Currently, returns a named list of the same length as pars where within each element of the list
#' is a numeric matrix (vector parameters) or vector (scalar parameters). For matrices, each row is a time point (or dimension) of the parameter and each column
#' is a named summary. The names are accessible with `colnames`. For vectors (scalar parameters), each element is a named summary.
#'
#'
#' @examples
#'
#' #' signal = c(rep(0, 50), rep(10, 50))
#' noise = rep(1, 100)
#' noise_var = rep(1, 100)
#' for (k in 2:100){
#'   noise_var[k] = exp(0.9*log(noise_var[k-1]) + rnorm(1, 0, 0.5))
#'   noise[k] = rnorm(1, 0, sqrt(noise_var[k])) }
#'
#' y = signal + noise
#' model_spec = dsp_spec(family = "gaussian", model = "changepoint")
#' mcmc_output = dsp_fit(y, model_spec = model_spec)
#' summary(mcmc_output)
#'
#' @export

summary.dsp <- function(object, pars, probs = c(0.025, 0.25, 0.50, 0.75, 0.975), ...){

  # supplied or default pars?

  if(missing(pars)){
    pars <- names(object$mcmc_output)
  }else{
    pars_in <- pars %in% names(object$mcmc_output)

    if(all(!pars_in)) stop("None of the requested parameters are available in object$mcmc_output.")

    if(!all(pars_in)){

      warning(paste("The following parameters were requested, but not available in model output:", paste(pars[!pars_in], collapse = ", ")))

      pars <- pars[pars_in]

    }

  }

  col_names <- c("mean", paste0("q_", round(100*probs, 1)))
  # Calculate the mean and quantiles, but return a list of same length, etc

  out_list <- purrr::map(object$mcmc_output[pars], \(samps){

    if(!is.null(dim(samps))){

      apply(samps, MARGIN = 2, summary_fun) |>
        t()

    }else{

      summary_fun(samps)

    }

  })

  return(out_list)

}

summary_fun <- function(col) {
  nums <- c(mean(col), quantile(col, probs = probs))
  names(nums) <- col_names
  nums
}

#----------------------------------------------------------------------------
#' Print a summary of dsp object
#'
#' @inheritParams summary.dsp
#'
#' @details
#' A brief summary of the settings used to fit the model including number of iterations,
#' burn in, and thinning rates.
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
#' print(mcmc_output)
#'
#'
#'
#' @export

print.dsp <- function(object, ...){

  # TODO add print for model spec here?
  cat("Total number of MCMC samples burned in:", object$mcpar["nburn"], "\n")
  cat("Thinning interval use:", object$mcpar["nskip"], "\n")
  cat("Total number of MCMC samples saved:", object$mcpar["nsave"])


  invisible(NULL)

}


