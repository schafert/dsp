#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import rlang
#' @importFrom fda create.bspline.basis
#' @importFrom fda eval.basis
#' @importFrom glue glue
#' @importFrom graphics legend
#' @importFrom lifecycle deprecated
#' @importFrom MCMCpack rinvgamma
#' @importFrom mgcv rig
#' @importFrom purrr map
#' @importFrom RcppZiggurat zrnorm
#' @importFrom spam chol
#' @importFrom stats dgamma
#' @importFrom stats dnbinom
#' @importFrom stats dpois
#' @importFrom stats rnbinom
#' @importFrom truncdist rtrunc
## usethis namespace: end
NULL

# Internal package state. Counts the state draws in which the precision matrix
# was too ill-conditioned to factorize and sampleBTF() fell back to flooring the
# evolution variances. dsp_fit() resets this before a fit and reports the total
# afterwards.
.dsp_state <- new.env(parent = emptyenv())
.dsp_state$n_illcond <- 0L

#' Dynamic Shrinkage Process (dsp) Object
#'
#' Check whether an object is a fitted \code{dsp} model.
#'
#' @details
#' A \code{dsp} object is a clcode{\link{dsp_fit}}.
#' It contains \code{mcmc_output}, a list of posterior draws whose matrix
#' elements are \code{nsave} br elements are of
#' length \code{nsave}; \code{DIC}, the deviance information criterion and
#' effective number of parametber of saved draws,
#' the burn-in and the thinning interval; and \code{model_spec}, the
#' \code{dsp_spec} object desc
#'
#' @param object any \R object
#' @return Logical, whether \code{object} inherits from class \code{dsp}.
#' @aliases dsp-class
#' @export
is.dsp <- function(object){
  inherits(object, "dsp")
}
