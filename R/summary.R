#----------------------------------------------------------------------------
#' Summarize a fitted dsp model
#'
#' Produces a compact summary of a fitted model: the model specification, the
#' MCMC settings, the deviance information criterion and one line per
#' parameter.
#'
#' @param object object of class dsp from [dsp_fit()]
#' @param pars parameter names to summarize; defaults to all parameters
#'   returned by [par_names()]
#' @param probs numeric vector of length two giving the probabilities of the
#'   posterior interval reported for scalar parameters. Defaults to
#'   c(0.025, 0.975)
#' @param ess logical; if TRUE, also report the minimum and median pointwise
#'   effective sample size of each time-varying or vector parameter. Defaults
#'   to FALSE, because this requires one spectral estimate per time point and
#'   takes several seconds for long series.
#' @param ... currently not used
#'
#' @details
#' Scalar parameters are summarized by their posterior mean, standard
#' deviation, interval and effective sample size. A parameter indexed by time
#' (or by basis function, or by predictor) is summarized in one line by its
#' length and the range of its pointwise posterior means, and, if
#' \code{ess = TRUE}, the minimum and median of its pointwise effective sample
#' sizes. A regression coefficient array \code{beta} is reported one line per
#' predictor.
#'
#' Effective sample sizes are computed by [coda::effectiveSize()]. Elements
#' with no posterior variation, such as the fixed initial values of the
#' evolution variances, are excluded from the minimum and median.
#'
#' Pointwise posterior summaries of a time-varying parameter are returned by
#' [fitted.dsp()], and the posterior draws by [as.mcmc.dsp()], for which the
#' diagnostic functions in \pkg{coda} are available.
#'
#' @returns
#' An object of class \code{summary.dsp}, a list with components
#' \describe{
#'   \item{model_spec}{the [dsp_spec()] object of the fitted model}
#'   \item{mcpar}{the number of saved draws, the burn-in and \code{nskip}}
#'   \item{nobs}{the number of time points}
#'   \item{DIC}{the data frame returned by [DIC()], or \code{NULL} if the DIC
#'     was not computed}
#'   \item{scalars}{a matrix with one row per scalar parameter and columns
#'     \code{mean}, \code{sd}, the two quantiles and \code{ess}, or
#'     \code{NULL} if there are none}
#'   \item{vectors}{a data frame with one row per non-scalar parameter and
#'     columns \code{length}, \code{mean_min} and \code{mean_max}, and, if
#'     \code{ess = TRUE}, \code{ess_min} and \code{ess_median}; or \code{NULL}
#'     if there are none}
#' }
#'
#' @seealso [fitted.dsp()], [DIC()], [as.mcmc.dsp()]
#'
#' @examples
#' set.seed(200)
#' signal = c(rep(0, 50), rep(10, 50))
#' noise = rep(1, 100)
#' noise_var = rep(1, 100)
#' for (k in 2:100){
#'   noise_var[k] = exp(0.9*log(noise_var[k-1]) + rnorm(1, 0, 0.5))
#'   noise[k] = rnorm(1, 0, sqrt(noise_var[k])) }
#'
#' y = signal + noise
#' model_spec = dsp_spec(family = "gaussian", model = "changepoint",
#'                       D = 1, useAnom = TRUE, obsSV = "SV")
#' mcmc_output = dsp_fit(y, model_spec = model_spec, nsave = 200, nburn = 200)
#'
#' summary(mcmc_output)
#'
#' # pointwise summaries of the mean function
#' head(fitted(mcmc_output))
#'
#' @method summary dsp
#' @export

summary.dsp <- function(object, pars, probs = c(0.025, 0.975), ess = FALSE, ...){

  check_dsp(object)

  if(!is.numeric(probs) || length(probs) != 2 || anyNA(probs) ||
     any(probs < 0 | probs > 1))
    stop("'probs' must be two probabilities between 0 and 1.", call. = FALSE)

  draws <- get_draws(object)

  # supplied or default pars?

  if(missing(pars)){
    pars <- names(draws)
  }else{
    pars_in <- pars %in% names(draws)

    if(all(!pars_in)) stop("None of the requested parameters are available in the fitted model.", call. = FALSE)

    if(!all(pars_in)){

      warning(paste("The following parameters were requested, but not available in model output:", paste(pars[!pars_in], collapse = ", ")))

      pars <- pars[pars_in]

    }

  }

  # Split each parameter into blocks of draws. Every block becomes one line of
  # the summary: a scalar parameter is one block, a vector parameter is one
  # block, and an nsave x T x J array is J blocks, one per predictor.
  blocks <- list()
  for(nm in pars){
    samps <- draws[[nm]]
    d <- dim(samps)
    if(length(d) == 3){
      for(j in seq_len(d[3]))
        blocks[[paste0(nm, "[, ", j, "]")]] <- matrix(samps[, , j], nrow = d[1])
    }else{
      blocks[[nm]] <- matrix(as.vector(samps), nrow = NROW(samps))
    }
  }

  is_scalar <- vapply(blocks, \(b) ncol(b) == 1, logical(1))

  scalars <- do.call(rbind, lapply(names(blocks)[is_scalar], \(nm){
    b <- blocks[[nm]]
    c(mean = mean(b), sd = stats::sd(b), stats::quantile(b, probs),
      ess = ess_varying(b))
  }))
  if(!is.null(scalars)) rownames(scalars) <- names(blocks)[is_scalar]

  vectors <- do.call(rbind, lapply(names(blocks)[!is_scalar], \(nm){
    means <- colMeans(blocks[[nm]])
    out <- data.frame(length = length(means),
                      mean_min = min(means), mean_max = max(means),
                      row.names = nm)
    if(ess){
      e <- ess_varying(blocks[[nm]])
      out$ess_min <- suppressWarnings(min(e, na.rm = TRUE))
      out$ess_median <- stats::median(e, na.rm = TRUE)
    }
    out
  }))

  structure(list(model_spec = object$model_spec,
                 mcpar = object$mcpar,
                 nobs = nobs(object),
                 DIC = if(is.null(object$DIC$DIC)) NULL else DIC(object),
                 scalars = scalars,
                 vectors = vectors),
            class = "summary.dsp")

}

#----------------------------------------------------------------------------
#' Effective sample sizes of the varying columns of a draws matrix
#'
#' @param b numeric matrix of posterior draws, one column per element
#' @return A numeric vector of length \code{ncol(b)}, \code{NA} for a column
#'   with no posterior variation, for which [coda::effectiveSize()] reports 0.
#' @keywords internal

ess_varying <- function(b){

  varies <- apply(b, 2, stats::var) > 0
  out <- rep(NA_real_, ncol(b))
  if(any(varies)) out[varies] <- coda::effectiveSize(b[, varies, drop = FALSE])
  unname(out)

}

#----------------------------------------------------------------------------
#' Print a summary of a fitted dsp model
#'
#' @param x object of class \code{summary.dsp} from [summary.dsp()]
#' @param digits number of significant digits to print
#' @param ... currently not used
#' @returns \code{x}, invisibly.
#' @method print summary.dsp
#' @export

print.summary.dsp <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  # Format each number on its own, so that one large value does not push a
  # whole column into scientific notation
  fmt <- \(z) vapply(z, \(zi) format(signif(zi, digits)), character(1))

  print(x$model_spec)

  # nskip iterations are discarded between saved draws, so the thinning
  # interval, as coda::thin() reports it, is nskip + 1
  cat("\nMCMC:", x$mcpar["nsave"], "saved draws (burn-in",
      paste0(x$mcpar["nburn"], ", thinning interval ", x$mcpar["nskip"] + 1, ");"),
      "T =", x$nobs, "\n")

  if(!is.null(x$vectors)){
    cat("\nTime-varying and vector parameters:\n")
    v <- x$vectors
    tab <- data.frame(length = v$length,
                      `posterior mean range` = paste0("[", fmt(v$mean_min), ", ", fmt(v$mean_max), "]"),
                      row.names = rownames(v), check.names = FALSE)
    if(!is.null(v$ess_min)){
      tab$`ESS min` <- round(v$ess_min)
      tab$`ESS median` <- round(v$ess_median)
    }
    print(tab, right = TRUE)
  }

  if(!is.null(x$scalars)){
    cat("\nScalar parameters:\n")
    s <- x$scalars
    tab <- apply(s[, colnames(s) != "ess", drop = FALSE], 2, fmt)
    tab <- matrix(tab, nrow = nrow(s),
                  dimnames = list(rownames(s), setdiff(colnames(s), "ess")))
    tab <- cbind(tab, ess = round(s[, "ess"]))
    print(noquote(tab), right = TRUE)
  }

  if(!is.null(x$DIC)){
    # DIC is compared across fits by differences, so print it to a fixed
    # number of decimals rather than significant digits
    fmt_dic <- \(z) format(round(z, 1), nsmall = 1)
    cat("\nDIC:", fmt_dic(x$DIC$DIC[1]), paste0("(p_D = ", fmt_dic(x$DIC$p_eff[1]), ");"),
        "with p_V:", fmt_dic(x$DIC$DIC[2]), paste0("(p_V = ", fmt_dic(x$DIC$p_eff[2]), ")\n"))
  }

  invisible(x)

}
