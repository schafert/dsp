#----------------------------------------------------------------------------
#' Detect changepoints from a changepoint model
#'
#' Computes, for each time point, the posterior probability that the
#' increment of the mean exceeds the threshold of the changepoint model, and
#' reports the time points at which this probability is at least `cp_thres`.
#'
#' @param object object of class dsp from [dsp_fit()], fitted with
#'   `model = "changepoint"`.
#' @param cp_thres (default 0.5) cutoff for the posterior probability of
#'   threshold exceedance needed to label a time point a changepoint.
#' @param times optional vector of length `T` giving the observation times,
#'   such as dates. If supplied, a `time` column is added to the result.
#' @param all logical; if `FALSE` (default), return only the detected
#'   changepoints. If `TRUE`, return one row for every time point.
#'
#' @details
#' The changepoint model of Wu et al. (2024) uses a thresholding mechanism with
#' a latent indicator variable. A time point \eqn{t} is a candidate
#' changepoint when the squared increment \eqn{\omega_t^2} exceeds the
#' threshold \eqn{\exp(\gamma)}. The posterior probability of this event is
#' estimated by the proportion of saved draws in which it occurs. Because
#' \eqn{\omega} is the `D`th difference of \eqn{\mu}, no increment is defined
#' at the first `D` time points, and their probability is `NA`.
#'
#' @returns
#' A data frame with columns `index` (the position in the observed series),
#' `time` (only if `times` is supplied), and `probability` (the posterior
#' probability of threshold exceedance). If `all = FALSE`, it has one row per
#' detected changepoint, and zero rows if none is detected. If `all = TRUE`,
#' it has one row per time point and an additional logical column
#' `changepoint`.
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
#'                       D = 1, useAnom = TRUE)
#' mcmc_output = dsp_fit(y, model_spec = model_spec, nsave = 200, nburn = 200)
#' changepoints(mcmc_output)
#' changepoints(mcmc_output, cp_thres = 0.8,
#'              times = seq(as.Date("2000-01-01"), by = "month", length.out = 100))
#' head(changepoints(mcmc_output, all = TRUE))
#'
#' @seealso [plot.dsp()], which marks the detected changepoints, and
#'   [credible_bands()], which assesses the increments jointly.
#' @export
changepoints <- function(object, cp_thres = 0.5, times = NULL, all = FALSE){

  check_dsp(object)

  if(object$model_spec$model != "changepoint")
    stop("changepoints() requires a model fitted with model = \"changepoint\".",
         call. = FALSE)

  if(!is.logical(all) || length(all) != 1L || is.na(all))
    stop("all must be TRUE or FALSE.", call. = FALSE)

  probability <- cp_probability(object, cp_thres)
  n <- length(probability)

  if(!is.null(times) && length(times) != n)
    stop(sprintf("times must have length %d, the length of the series.", n),
         call. = FALSE)

  detected <- !is.na(probability) & probability >= cp_thres
  rows <- if(all) seq_len(n) else which(detected)

  out <- data.frame(index = as.integer(rows))
  if(!is.null(times)) out$time <- times[rows]
  out$probability <- probability[rows]
  if(all) out$changepoint <- detected
  out
}

# Posterior probability of threshold exceedance at each of the T time points,
# with NA at the first D, where no increment is defined.
cp_probability <- function(object, cp_thres){

  if(!all(c("omega", "gamma") %in% names(object$mcmc_output)))
    stop("Can't calculate changepoint evidence without samples of omega and gamma.",
         call. = FALSE)

  if(!is.numeric(cp_thres) || length(cp_thres) != 1L || is.na(cp_thres) ||
     cp_thres < 0 || cp_thres > 1)
    stop("cp_thres must be a numeric value between 0 and 1.", call. = FALSE)

  D <- object$model_spec$arguments$D
  cp_mat <- apply(object$mcmc_output$omega^2, MARGIN = 2,
                  \(samp){samp > exp(object$mcmc_output$gamma)})
  c(rep(NA_real_, D), colMeans(cp_mat))
}
