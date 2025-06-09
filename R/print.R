#' Print a summary of a model specification or fitted dsp object
#'
#' Prints basic information about the MCMC sampling settings for the fitted model
#'
#' @param object object of class dsp from [dsp_fit()]
#' @param ... currently not used
#'
#' @details
#' A brief summary of the settings used to fit the model including number of iterations,
#' burn in, and thinning rates.
#'
#' @returns NULL
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
#' @method print dsp
#' @export
#'

print.dsp <- function(object, ...){

  print(object$model_spec)
  cat("\nTotal number of MCMC samples burned in:", object$mcpar["nburn"], "\n")
  cat("Thinning interval use:", object$mcpar["nskip"], "\n")
  cat("Total number of MCMC samples saved:", object$mcpar["nsave"])


  invisible(NULL)

}


#' Print a summary of a model specification or fitted dsp object
#'
#' Prints basic information about the model specification
#'
#' @param object object of class dsp_spec from [dsp_spec()]
#' @param ... currently not used
#'
#' @returns NULL
#' @export
#'
#' @examples
#' model_spec = dsp_spec(family = "gaussian", model = "changepoint")
#' print(model_spec)
#'
#' @method print dsp_spec
#' @rdname print.dsp
print.dsp_spec <- function(object, ...){

  if(object$family == "gaussian"){

    cat("Gaussian likelihood with ")

    switch(object$model,
           changepoint = cat("changepoint and outlier detection, "),
           smoothing = cat("Bayesian trend smoothing, "),
           regression = cat("time-varying regression, "),
           bspline = cat("B-spline smoothing splines, "))

    switch(object$arguments$obsSV,
           const = cat("constant error variance, and "),
           SV = cat("stochastic volatility, and "),
           ASV = cat("adaptive stochastic volatility, and "))

    cat(object$arguments$D, "degree of differencing.")
  }

  if(object$family == "negbinom"){
    cat("Negative binomial likelihood with", object$arguments$D, "degree of differencing.")
  }

  invisible(NULL)

}
