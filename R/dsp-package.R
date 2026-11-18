#' @title dsp: Dynamic shrinkage process and change point detection
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

# Internal count of successful square-root recoveries across the six Gaussian
# state/coefficient samplers. dsp_fit() resets it before each fit. Recovery
# never floors variances or repairs the precision matrix.
.dsp_state <- new.env(parent = emptyenv())
.dsp_state$n_sqrt <- 0L

#' Dynamic shrinkage process (dsp) object
#'
#' Check whether an object is a fitted \code{dsp} model.
#'
#' @details
#' A \code{dsp} object is a classified list created by \code{\link{dsp_fit}}.
#' It contains \code{mcmc_output}, a list of posterior draws whose matrix
#' elements are \code{nsave} by \code{T} and whose vector elements are of
#' length \code{nsave}; \code{DIC}, the deviance information criterion and
#' effective number of parameters; \code{mcpar}, the number of saved draws,
#' the burn-in and the thinning interval; and \code{model_spec}, the
#' \code{dsp_spec} object describing the fitted model.
#' The \code{y} component contains the original observations, including missing values.
#'
#' The accessor functions [par_names()], [fitted.dsp()], [DIC()],
#' [get_spec()] and [nobs.dsp()] extract these components without referring
#' to the structure of the object.
#'
#' \code{coef()} returns the posterior means of the regression coefficients for
#' \code{model = "regression"} or the basis coefficients for \code{model = "bspline"}.
#' \code{residuals()} subtracts the posterior mean of \eqn{\mu_t} from \eqn{y_t}.
#' \code{confint(object, pars = "mu", level = 0.95)} returns pointwise credible
#' intervals for the specified parameter.
#' \code{logLik()} returns the posterior mean log-likelihood. Set
#' \code{summary = "draws"} to return the value at each saved iteration.
#' \code{predict()} returns the same summaries as \code{fitted(object, pars = "mu")}
#' for regression, smoothing, and B-spline models. For changepoint models, it
#' returns detected changepoint indices; \code{cp_prop = TRUE} also returns
#' their probabilities.
#'
#' @seealso [dsp-accessors], [as.mcmc.dsp()], [summary.dsp()]
#'
#' @param object any \R object
#' @return Logical, whether \code{object} inherits from class \code{dsp}.
#' @aliases dsp-class
#' @export
is.dsp <- function(object){
  inherits(object, "dsp")
}
