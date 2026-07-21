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
