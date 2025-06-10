#' Model Specification
#'
#' Method for creating dsp specification object prior to fitting.
#'
#' @param family A character string specifying the model family. Must be one of:
#'   \itemize{
#'     \item "gaussian": Gaussian family.
#'     \item "negbinom": Negative binomial family.
#'   }
#' @param model A character string specifying the model type:
#'   \itemize{
#'     \item \code{family} = "gaussian":
#'       \itemize{
#'         \item "changepoint": Change point detection with
#'         Adaptive Bayesian Changepoint analysis and local Outlier (ABCO),
#'         \item "smoothing": Bayesian smoothing,
#'         \item "regression": Time-varying regression,
#'         \item "bspline": Bayesian smoothing with B-spline for irregularly spaced or functional time-series.
#'       }
#'     \item \code{family} = "negbinom":
#'       \itemize{
#'         \item "smoothing": Bayesian smoothing.
#'       }
#'   }
#' @param ... Additional arguments based on \code{family} and \code{model}.
#' For \code{family} = "gaussian" and \code{model} = "regression", \code{X} is required.
#' Other arguments are optional and default to predefined value.
#'
#' \tabular{lll}{
#'  **Family**    \tab  **Model**       \tab    **Required Arguments**  \cr
#'  "gaussian"    \tab "changepoint"     \tab \code{D}, \code{useAnom}, \code{obsSV} \cr
#'  "gaussian"    \tab "smoothing"       \tab \code{D}, \code{evol_error}, \code{obsSV}, \code{zero_error}, \code{D_asv}, \code{evol_error_asv}, \code{nugget_asv}\cr
#'  "gaussian"    \tab "regression"      \tab \code{D}, \code{evol_error}, \code{obsSV}, \code{X}, \code{D_asv}, \code{evol_error_asv}, \code{nugget_asv} \cr
#'  "gaussian"    \tab "bspline"        \tab \code{D}, \code{evol_error}, \code{obsSV}, \code{times}, \code{num_knots}  \cr
#'  "negbinom"    \tab "smoothing"       \tab \code{D}, \code{evol_error}, \code{r_init}, \code{r_sample}, \code{offset}  \cr
#' }
#'
#'   \itemize{
#'     \item **Shared Arguments**:
#'       \itemize{
#'         \item \code{D}: integer scalar indicating degree of differencing.
#'         Implementation is available D = 0, 1, or 2 (default) for \code{family} = "gaussian",
#'         and D = 1 and 2 (default) for \code{family} = "negbinomial".
#'         \item \code{obsSV}:Options for modeling the error variance for \code{family} = "gaussian".
#'         It must be one of the following:
#'            \itemize{
#'              \item "const": Constant error variance for all time points (default);
#'              \item "SV": Stochastic Volatility model;
#'              \item "ASV": Adaptive Stochastic Volatility model;
#'            }
#'         \item \code{D_asv}: integer; degree of differencing (0, 1 (default), or 2) for the ASV model. Only used when \code{obsSV = "ASV"}.
#'         \item \code{evol_error_asv}: character; "HS" by default. evolution error distribution for the ASV model. Must be one of the five options used in \code{evol_error}. Only used when \code{obsSV = "ASV"}
#'         \item \code{nugget_asv}: logical; if \code{TRUE} (default), fits the nugget variant of the ASV model. Only used when \code{obsSV = "ASV"}.
#'         \item \code{evol_error}: the evolution error distribution; .
#'            \itemize{
#'              \item the dynamic horseshoe prior ("DHS") (default);
#'              \item the static horseshoe prior ("HS");
#'              \item the Bayesian lasso ("BL") (\code{family} = "gaussian" only);
#'              \item the normal stochastic volatility model ("SV") (\code{family} = "gaussian" only);
#'              \item the normal-inverse-gamma prior ("NIG") (\code{family} = "gaussian" only).
#'            }
#'       }
#'     \item **Model-specific Arguments**:
#'      \itemize{
#'        \item \code{model} = "changepoint":
#'          \itemize{
#'            \item \code{useAnom}: logical; Defaults to FALSE. if TRUE, include an anomaly component
#'            in the observation equation.
#'          }
#'        \item \code{family} = "gaussian", \code{model} = "smoothing":
#'        \itemize{
#'          \item \code{zero_error}: one of the 5 distributions listed in \code{"evol_error"};
#'          Additional penalty distribution for the conditional expectation i.e., shrinkage to zero (optional).
#'        }
#'        \item \code{family} = "gaussian", \code{model} = "regression":
#'        \itemize{
#'          \item \code{X}: matrix or data frame; T x p matrix of time series predictors
#'        }
#'        \item \code{model} = "bspline":
#'        \itemize{
#'          \item \code{times}: length \code{T} vector; observation indices; if NULL (Default), assume equally spaced.
#'          \item \code{num_knots}: numeric; the number of knots to be used for bspline. defaults to 20.
#'        }
#'        \item \code{family} = "negbinomial", \code{model} = "smoothing":
#'        \itemize{
#'          \item \code{r_init}: numeric; initial value (defaults to 5) of MCMC sampling for
#'          overdispersion parameter.
#'          \item \code{r_sample}: logical; If `TRUE`, the overdispersion is sampled,
#'          if `FALSE` (default), fixed at "\code{r_init}.
#'          \item \code{offset}: length \code{T} vector; offset values for the log conditional expectation
#'        }
#'     }
#'   }
#'
#' @return A list containing the model specification.
#'
#' @examples
#' model_spec <- dsp_spec(family = "gaussian",
#'                       model = "changepoint")
#'
#'
#' @export
dsp_spec <- function(family,
                     model,
                     ...){
  if(!family %in% c("gaussian", "negbinomial")) stop('family must be gaussian or negbinomial')
  input_args <- list(...)
  if(family == "gaussian"){
    if(!model %in% c("changepoint", "regression",
                     "smoothing","bspline")){
      stop('model must be "changepoint","regression","smoothing" or "bspline for gaussian family')
    }
  }else{
    if(model != "smoothing"){
      stop('model must be "smoothing" for negative binomial family')
    }
  }
  # Extract obsSV if provided
  obsSV <- input_args$obsSV
  if (is.null(obsSV)) obsSV <- "NONE"  # default fallback if needed

  #inputted arguments
  required_args <- function(family,model){
    if(family == "gaussian"){
      if(model == "changepoint"){
        return(c("D", "useAnom", "obsSV"))
      }else if(model == "smoothing" && obsSV == "ASV"){
        return(c("evol_error","D","obsSV","zero_error","D_asv","evol_error_asv","nugget_asv"))
      }else if(model == "smoothing" && obsSV != "ASV"){
        return(c("evol_error","D","obsSV","zero_error"))
      }else if(model == "regression" && obsSV == "ASV"){
        return(c("X", "evol_error", "D","obsSV","D_asv","evol_error_asv","nugget_asv"))
      }else if(model == "regression" && obsSV == "ASV"){
        return(c("X", "evol_error", "D","obsSV"))
      }else{
        return(c("evol_error", "D","times","num_knots"))
      }
    }else{
      # negative binomial smoothing
      return(c("evol_error","D","r_init","r_sample","offset"))
    }
  }
  # Handling unnecessary arguments supplied as inputs
  extraArgs = setdiff(names(input_args),required_args(family,model))
  if(length(extraArgs) > 0){
    warning(
      sprintf("The following arguments are not used for %s family %s model and will be ignored: %s",
              family,model, paste(extraArgs, collapse = ", ")),
      call. = FALSE)
    for(arg in extraArgs){
      input_args[[arg]] = NULL
    }
  }
  # check the input type and throw error if the inputs are invalid
  validation_rules <- list(
    X = function(family,x) { if (any(is.na(x)) || !(is.matrix(x) || is.data.frame(x))) stop("X must be a dataframe or a matrix")},
    D = function(family,x) {
        if(family == "gaussian"){
          if (!x %in% c(0, 1, 2)) stop("D must be 0, 1 or 2 for gaussian models.")
        }else{
          if (!x %in% c(1, 2)) stop("D must be 1 or 2 for negative binomial models")
        }
      },
    D_asv = function(family,x) {
      if (!x %in% c(0, 1, 2)) stop("D must be 0, 1 or 2 for gaussian models.")
    },
    nugget_asv = function(family,x) {
      if (is.na(x)|| !is.logical(x)) stop("useAnom must be TRUE or FALSE.")
    },
    evol_error_asv = function(family,x) {
      if (!x %in% c("HS", "DHS", "NIG", "SV", "BL")) {
        stop("evol_error_asv must be one of 'HS', 'DHS', 'NIG', 'SV', 'BL' for gaussian models")
      }
    },
    evol_error = function(family,x) {
      if(family == "gaussian"){
        if (!x %in% c("HS", "DHS", "NIG", "SV", "BL")) {
          stop("evol_error must be one of 'HS', 'DHS', 'NIG', 'SV', 'BL' for gaussian models")
        }else{
          if (!x %in% c("HS", "DHS", "NIG", "SV", "BL")) {
            stop("evol_error must be one of 'HS', 'DHS' for negative binomial models" )
          }
        }
      }
    },
    useAnom = function(family,x) { if (is.na(x)|| !is.logical(x)) stop("useAnom must be TRUE or FALSE.") },
    obsSV = function(family,x) { if (!x %in% c("const", "SV", "ASV")) stop("obsSV must be one of 'const', 'SV', 'ASV'.") },
    zero_error = function(family,x) { if (!x %in% c("HS", "DHS", "NIG", "SV", "BL")) stop("zero_error must be one of 'HS', 'DHS', 'NIG', 'SV', 'BL'.")},
    times = function(family,x) { if (any(is.na(x)) || !is.numeric(x) || any(x < 0 || x != as.integer(x))) stop("must be a vector of positive integer or NULL") },
    num_knots = function(family,x) { if (is.na(x) || !is.numeric(x) || x<0 || x != as.integer(x))  stop("num_knots must be a positive integer.") },
    r_init = function(family,x){ if (!is.numeric(x) || x < 0) stop("r_init nust be a positive number.")},
    r_sample = function(family,x) { if (is.na(x) || !is.logical(x)) stop("r_sample must be TRUE or FALSE.")},
    offset = function(family,x){ if (any(is.na(x)) || !is.numeric(x)) stop("offset must be a vector") }
  )
  if(length(input_args) >0){
    for(arg in names(input_args)){
      tryCatch(
        validation_rules[[arg]](family,input_args[[arg]]),
        error = function(e) { stop(sprintf("invalid input for %s model in %s family: %s",
                                           model, family, e$message)) }
      )
    }
  }
  # If any of the required arguments are missing give them default values
  default_args <- list(D =2, useAnom = FALSE, obsSV = "const",
                       evol_error = "DHS", zero_error = NULL, num_knots = 20,
                       r_init = 5, r_sample = FALSE, offset = 0,
                       D_asv = 1,evol_error_asv = "HS",nugget_asv = TRUE)
  requiredArgs = setdiff(required_args(family,model), names(input_args))
  # X is required to run regression
  if("X" %in% requiredArgs) stop("To run regression, the design matrix 'X' needs to be specified.")
  if(length(requiredArgs)>0){
    for(arg in requiredArgs){
      input_args[[arg]] = default_args[[arg]]
      message(sprintf("Argument '%s' is missing. Using default value: %s", arg, default_args[[arg]]))
    }
  }
  ret = structure(list(
    family = family,
    model = model,
    arguments = input_args
  ),
  class = c("dsp_spec"))
  return(ret)

}

#' MCMC Sampler for Models with Dynamic Shrinkage Processes
#'
#' Wrapper function for fitting models with Dynamic Shrinkage Processes (DSP), including:
#' \itemize{
#'         \item Adaptive Bayesian Changepoint analysis and local Outlier (ABCO),
#'         \item Bayesian Trend Filter for Gaussian Data
#'         \item Time-varying Regression
#'         \item Bayesian Trend Filter with B-spline for irregularly spaced or functional time-series.
#'         \item Bayesian Smoothing for Count Data
#'}
#'
#' @param y a numeric vector of the \code{T x 1} vector of time series observations
#' @param model_spec a list containing model specification generated from [dsp_spec()].
#' @param nsave integer scalar (default = 1000); number of MCMC iterations to record
#' @param nburn integer scalar (default = 1000); number of MCMC iterations to discard (burn-in)
#' @param nskip integer scalar (default = 4); number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param computeDIC logical; if TRUE (default), compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should extra information on progress be printed to the console? Defaults to FALSE
#' @param ... optional additional arguments to pass to the MCMC sampler.
#'
#' @return \code{dsp_fit} returns an object of class "\code{dsp}".
#'
#' An object of class "\code{dsp}" is defined as a list containing at least the following components:
#'    \item{mcmc_output}{a list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}}
#'    \item{DIC}{Deviance Information Criterion}
#'    \item{mcpar}{named vector of supplied nsave, nburn, and nskip}
#'    \item{model_spec}{the object supplied for model_spec argument}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues when family is "gaussian".
#'
#' @examples
#' \dontrun{
#' beta <- make_signal(name = "bumps", n = 300)
#' y <- rnbinom(n = length(beta), size = 5, mu = beta)
#' # Creating a `model_spec` object
#' model_spec <- dsp_spec(
#'   family = "negbinomial",
#'   model = "smoothing")
#'
#' fit <- dsp_fit(
#'   y = y,
#'   model_spec = model_spec,
#'   nburn = 5000,
#'   nsave = 5000)
#'}
#'
#' @export
dsp_fit = function(y, model_spec,
                      nsave = 1000, nburn = 1000, nskip = 4,
                      computeDIC = TRUE,
                      verbose = TRUE,
                      ...){
  if (!inherits(model_spec, "dsp_spec")) {
    stop("'model_spec' must be an object created by 'dsp_spec()'.")
  }

  family = model_spec$family
  model = model_spec$model
  if(!is.null(model_spec$arguments$zero_error)){
    model = "smoothing_sparse"
  }
  fitter_list <- list(
    gaussian = list(
      changepoint = abco, #changepoint model,
      regression = btf_reg, #Bayesian Trend Filter with Regression,
      smoothing = btf, #Bayesian Trend Filter,
      smoothing_sparse = btf_sparse, #Bayesian Trend Filter (more shrinkage),
      bspline = btf_bspline #Bayesian Trend Filter with Bspline
    ),
    negbinomial = list(
      smoothing = btf_nb #  Bayesian Trend Filter with Negative Binomial
    )
  )
  fitter = fitter_list[[family]][[model]]
  if (is.null(fitter)) {
    stop("No valid fitter function found for the specified family and model.")
  }
  ## This is a temporary case until we figure out what to do about the DIC thing

  if(model == "changepoint"){
    input_args <- c(model_spec$arguments, list(y = y, nsave = nsave, nburn = nburn,
                                               nskip = nskip,
                                               verbose = verbose, ...))
  }else{
    input_args <- c(model_spec$arguments, list(y = y, nsave = nsave, nburn = nburn,
                                               nskip = nskip, computeDIC = computeDIC,
                                               verbose = verbose, ...))
  }

  mcmc_output = do.call(fitter, input_args)
  structure(list(mcmc_output = mcmc_output,
                 DIC = mcmc_output[c("DIC", "p_d")],
                 mcpar = c(nsave = nsave, nburn = nburn, nskip = nskip),
                 model_spec = model_spec),
            class = c(c("dsp"), mcmc_output$class))

}




