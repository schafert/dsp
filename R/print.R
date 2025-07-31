#' Print a summary of a model specification or fitted dsp object
#'
#' Method for printing basic information about the MCMC sampling settings for the fitted model
#'
#' @param x object of class dsp from [dsp_fit()]
#' @param ... currently not used
#'
#' @details
#' A brief summary of the settings used to fit the model including number of iterations,
#' burn in, and thinning rates.
#'
#' @returns NULL
#'
#' @examples
#' print(mcmc_output)
#'
#' @method print dsp
#' @rdname dsp_fit
#' @export
#'

print.dsp <- function(x, ...){

  print(x$model_spec)
  cat("\nTotal number of MCMC samples burned in:", x$mcpar["nburn"], "\n")
  cat("Thinning interval use:", x$mcpar["nskip"], "\n")
  cat("Total number of MCMC samples saved:", x$mcpar["nsave"])


  invisible(NULL)

}


#' Print a summary of a model specification or fitted dsp object
#'
#' Method for printing basic information about the model specification
#'
#' @param x object of class dsp_spec from [dsp_spec()]
#' @param ... currently not used
#'
#' @returns NULL
#'
#' @export
#'
#' @examples
#' print(model_spec)
#'
#' @method print dsp_spec
#' @rdname dsp_spec
print.dsp_spec <- function(x, ...){

  if(x$family == "gaussian"){

    cat("Gaussian likelihood with ")

    switch(x$model,
           changepoint = cat("changepoint and outlier detection, "),
           smoothing = cat("Bayesian trend smoothing, "),
           regression = cat("time-varying regression, "),
           bspline = cat("B-spline smoothing splines, "))

    switch(x$arguments$obsSV,
           const = cat("constant error variance, and "),
           SV = cat("stochastic volatility, and "),
           ASV = cat("adaptive stochastic volatility, and "))

    cat(x$arguments$D, "degree of differencing.")
  }

  if(x$family == "negbinom"){
    cat("Negative binomial likelihood with", x$arguments$D, "degree of differencing.")
  }

  invisible(NULL)

}
