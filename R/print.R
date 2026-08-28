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
  cat("Thinning interval used:", x$mcpar["nskip"], "\n")
  cat("Total number of MCMC samples saved:", x$mcpar["nsave"],"\n")


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

    if (x$model %in% c("changepoint", "smoothing", "regression")) {
      switch(x$arguments$obsSV,
        const = cat("constant error variance, and "),
        SV    = cat("stochastic volatility, and "),
        ASV   = cat("adaptive stochastic volatility, and "),
        cat("unspecified error variance model, and ")
      )
    }

    cat(x$arguments$D, if (isTRUE(x$arguments$D == 1)) "degree" else "degrees",
        "of differencing.\n")
  }

  if(x$family == "negbinomial"){
    cat("Negative binomial likelihood with", x$arguments$D,
        if (isTRUE(x$arguments$D == 1)) "degree" else "degrees", "of differencing. \n")
  }
  if (x$model == "changepoint"){
    cat("Prior for the mean function: Dynamic horseshoe prior with threshold autoregression\n")
  }
  if (!is.null(x$arguments$evol_error)) {
    cat("Prior for the mean function: ")
    switch(x$arguments$evol_error,
           HS  = cat("Horseshoe prior"),
           DHS = cat("Dynamic horseshoe prior"),
           NIG = cat("Normal-inverse-gamma prior"),
           BL  = cat("Bayesian lasso prior"),
           SV  = cat("Stochastic volatility prior"))
    cat("\n")
  }
  if(x$family == "gaussian"){
    cat("Prior for the variance function: ")
    if (x$model == "bspline"){
      cat("Normal-inverse-gamma prior\n")
    }else{
      switch(x$arguments$obsSV,
             const = cat("Normal-inverse-gamma prior"),
             SV = cat("Stochastic volatility prior"),
             ASV = {cat("Adaptive stochastic volatility with ")
               switch(x$arguments$evol_error_asv,
                      HS  = cat("Horseshoe prior"),
                      DHS = cat("Dynamic horseshoe prior"),
                      NIG = cat("Normal-inverse-gamma prior"),
                      BL  = cat("Bayesian lasso prior"),
                      SV  = cat("Stochastic volatility prior"))
             })
      cat("\n")
    }
  }

  invisible(NULL)
}
