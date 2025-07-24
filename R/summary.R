#----------------------------------------------------------------------------
#' Summarize DSP MCMC chains
#'
#' @param object object of class dsp from [dsp_fit()]
#' @param pars parameter names specified for summaries; currently defaults to all parameters named in object$mcmc_output
#' @param probs numeric vector of [quantile()]s requested for posterior summary of pars. Defaults to c(0.025, 0.25, 0.50, 0.75, 0.975)
#' @param ... currently not being used
#'
#' @returns
#' Returns a named list of the same length as pars where within each element of the list
#' is a numeric matrix (vector parameters) or vector (scalar parameters). For matrices, each row is a time point (or dimension) of the parameter and each column
#' is a named summary. The names are accessible with `colnames`. For vectors (scalar parameters), each element is a named summary.
#'
#'
#' @examples
#' set.seed(200)
#' mu <- make_signal(name = "quadratic", n = 300)        # Underlying trend with a ramp structure
#' var  <- make_signal(name = "bumps", n = 300)/2       # Time-varying standard deviation with bumps
#' y <- rnorm(n = 300, mean = mu, sd = sqrt(var))    # Observed data based on above.
#' # Specify DSP model with ASV observation noise and Horseshoe prior on evolution error
#' spec <- dsp_spec(family = "gaussian",
#'                  model = "smoothing",
#'                  obsSV = "ASV",
#'                  D_asv = 2)
#' # Fit the model (Note: longer MCMC runs may be required for stable inference)
#' fit <- dsp_fit(y, spec, nsave = 1000, nburn = 1000)
#' summary_fit <- summary(fit)
#' summary_fit$mu[,"mean"]  # contains estimated posterior mean
#' summary_fit$h[,"mean"] # conrtains log volatility
#'
#' @importFrom purrr map
#'
#' @method summary dsp
#' @export
summary.dsp <- function(object, pars, probs = c(0.025, 0.25, 0.50, 0.75, 0.975), ... ){

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

  out_list <- purrr::map(object$mcmc_output[pars], \(samps) {

    if (is.null(dim(samps))) {
      # Scalar
      return(summary_fun(samps, probs = probs, col_names = col_names))

    } else if (length(dim(samps)) == 2) {
      # Matrix (M x T)
      return(t(apply(samps, 2, summary_fun, probs = probs, col_names = col_names)))

    } else if (length(dim(samps)) == 3) {
      # Array (M x T x J)
      # Output: array of dim T x J x K (statistic)
      dims <- dim(samps)
      M <- dims[1]; T <- dims[2]; J <- dims[3]
      out <- array(NA, dim = c(T, J, length(col_names)),
                   dimnames = list(NULL, NULL, col_names))

      for (j in 1:J) {
        tmp <- apply(samps[,,j], 2, summary_fun, probs = probs, col_names = col_names)
        out[,j,] <- t(tmp)
      }
      return(out)

    } else {
      stop("Unsupported parameter shape: must be scalar, 2D, or 3D.")
    }
  })

  return(out_list)

}
#' @keywords internal
summary_fun <- function(col, probs, col_names) {
  nums <- c(mean(col), quantile(col, probs = probs))
  names(nums) <- col_names
  nums
}





