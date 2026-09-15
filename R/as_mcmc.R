#----------------------------------------------------------------------------
#' Coerce a fitted dsp model to an mcmc object
#'
#' Flattens the posterior draws stored in a \code{dsp} object into a single
#' matrix of class \code{\link[coda]{mcmc}}, so that the diagnostic and
#' plotting functions in \pkg{coda} can be applied to the fitted model.
#'
#' @param x an object of class `dsp`; see [dsp-class]
#' @param pars character vector of parameter names to retain; defaults to all
#'   parameters named in \code{x$mcmc_output}
#' @param ... Further arguments to be passed to specific methods
#'
#' @details
#' Each element of \code{x$mcmc_output} is stored with the posterior draws in
#' the first dimension, so that a scalar parameter is a vector of length
#' \code{nsave} and a time-varying parameter is an \code{nsave} by \code{T}
#' matrix or an \code{nsave} by \code{T} by \code{J} array. Every element is
#' flattened column by column and the results are bound into one matrix with
#' \code{nsave} rows. Columns are named by indexing the original parameter,
#' dropping the leading draw index. Column \code{"mu[3]"} therefore holds the
#' draws of the mean function at time 3, and column \code{"beta[3, 2]"} holds
#' the draws of the second regression coefficient at time 3.
#'
#' The deviance information criterion and the effective number of parameters
#' are posterior summaries rather than draws, and are dropped by the coercion.
#' They remain available in \code{x$DIC}.
#'
#' The sampler saves every \code{nskip + 1} iterations after burn-in, so the
#' returned object records a thinning interval of \code{nskip + 1}, a starting
#' iteration of \code{nburn + nskip + 1} and a final iteration of
#' \code{nburn + nsave * (nskip + 1)}.
#'
#'
#' @returns
#' An object of class \code{\link[coda]{mcmc}} with \code{nsave} rows and one
#' column per scalar element of the retained parameters named by \code{pars}.
#'
#' @seealso [dsp-class], [summary.dsp()], [coda::mcmc()]
#'
#' @examples
#' set.seed(200)
#' signal <- c(rep(0, 50), rep(10, 50))
#' y <- signal + rnorm(100)
#'
#' model_spec <- dsp_spec(family = "gaussian", model = "smoothing", D = 1)
#' fit <- dsp_fit(y, model_spec = model_spec, nsave = 200, nburn = 200)
#'
#' draws <- coda::as.mcmc(fit, pars = c("dhs_phi", "dhs_mean"))
#' summary(draws)
#' coda::effectiveSize(draws)
#'
#' @importFrom coda as.mcmc mcmc
#' @method as.mcmc dsp
#' @export

as.mcmc.dsp <- function(x, pars, ...){

  if(!is.dsp(x))
    stop("'x' must be an object of class 'dsp'.", call. = FALSE)

  nsave <- unname(x$mcpar["nsave"])
  nburn <- unname(x$mcpar["nburn"])
  nskip <- unname(x$mcpar["nskip"])

  # Retain only the elements that hold posterior draws. DIC and p_d are
  # posterior summaries, and any element whose leading dimension is not the
  # number of saved draws cannot be a monitored quantity.
  draws <- x$mcmc_output[!names(x$mcmc_output) %in% c("DIC", "p_d")]

  is_draw <- vapply(draws, function(samps){
    d <- dim(samps)
    if(is.null(d)) length(samps) == nsave else d[1] == nsave
  }, logical(1))

  draws <- draws[is_draw]

  if(length(draws) == 0)
    stop("No posterior draws were found in 'x$mcmc_output'.", call. = FALSE)

  # Supplied or default pars?
  if(missing(pars)){
    pars <- names(draws)
  } else {
    pars_in <- pars %in% names(draws)

    if(all(!pars_in))
      stop("None of the requested parameters are available in x$mcmc_output.",
           call. = FALSE)

    if(!all(pars_in))
      warning(paste("The following parameters were requested, but not available in model output:",
                    paste(pars[!pars_in], collapse = ", ")))

    pars <- pars[pars_in]
  }

  # Flatten each parameter to an nsave x prod(dim[-1]) matrix. as.vector() is
  # column-major, so the first index of the parameter varies fastest, which is
  # the order the column names are built in.
  flat <- lapply(pars, function(nm){
    samps <- draws[[nm]]
    out <- matrix(as.vector(samps), nrow = nsave)
    colnames(out) <- get_flatnames(nm, dim(samps)[-1])
    return(out)
  })

  out <- do.call(cbind, flat)

  thin <- nskip + 1

  return(coda::mcmc(out,
             start = nburn + thin,
             end   = nburn + nsave*thin,
             thin  = thin))

}

#----------------------------------------------------------------------------
#' Construct indexed column names for a flattened parameter
#'
#' Builds the column names used when a posterior array is flattened to a
#' matrix, in the column-major order produced by \code{as.vector()}.
#'
#' @param name character, the name of the parameter
#' @param dim integer vector of the parameter dimensions excluding the draws,
#'   or \code{NULL} for a scalar parameter
#' @param first_index integer, the value the indices start from
#' @return A character vector of length \code{prod(dim)}, or \code{name} itself
#'   when \code{dim} is \code{NULL} or empty.
#' @keywords internal

get_flatnames <- function(name, dim = NULL, first_index = 1){

  if(is.null(dim) || length(dim) == 0){
    return(name)
  }

  indices <- lapply(dim, \(ind){first_index:(ind + first_index - 1)})

  # expand.grid() varies the first index fastest, matching column-major order
  grid <- expand.grid(indices, KEEP.OUT.ATTRS = FALSE)

  return(paste0(name, "[", do.call(paste, c(grid, sep = ", ")), "]"))

}
