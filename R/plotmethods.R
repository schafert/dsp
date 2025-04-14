#----------------------------------------------------------------------------
#' Plot the Bayesian trend filtering fitted values
#'
#' Plot the BTF posterior mean of the conditional expectation with posterior credible intervals (pointwise and joint),
#' the observed data, and true curves (if known)
#'
#' @param mcmc_output an object of class "dsp" with parameter names 'mu' and/or 'yhat'
#' @param y the \code{T x 1} vector of time series observations
#' @param mu_true the \code{T x 1} vector of the true conditional mean
#' @param t01 the observation points; if NULL, assume \code{T} equally spaced points from 0 to 1
#' @param include_joint_bands logical; if TRUE, compute simultaneous credible bands
#'
#' @details
#'
#' The credible intervals plotted depends on whether `mcmc_output` contains predictions
#' of the data, i.e., it contains the name 'yhat'. If it does, then a plot is created with
#' credible intervals for the posterior prediction of y, else credible intervals are calculated
#' for the conditional expectation contained in the named slot 'mu'.
#'
#' @examples
#' signal = c(rep(0, 50), rep(10, 50))
#' noise = rep(1, 100)
#' noise_var = rep(1, 100)
#' for (k in 2:100){
#'   noise_var[k] = exp(0.9*log(noise_var[k-1]) + rnorm(1, 0, 0.5))
#'   noise[k] = rnorm(1, 0, sqrt(noise_var[k])) }
#'
#' y = signal + noise
#' model_spec = dsp_spec(family = "gaussian", model = "changepoint", mcmc_params = list('yhat', 'mu', "omega", "r"))
#' mcmc_output = dsp_fit(y, model_spec = model_spec)
#' plot(mcmc_output, y, mu_true = signal)
#'
#' @import coda
#' @export
plot.dsp = function(mcmc_output, y, mu_true = NULL, t01 = NULL, include_joint_bands = FALSE){

  # Time series:
  nT = length(y);
  if(is.null(t01)) t01 = seq(0, 1, length.out=nT)

  if(!("mu" %in% names(mcmc_output) || "yhat" %in% names(mcmc_output))){

    stop("Expected entries named 'mu' and 'yhat' in mcmc_output")

  }

  mu = colMeans(mcmc_output$mu)
  postY = mcmc_output$yhat


  # Credible intervals/bands:
  # TODO: restructure so that there's more conditioning because legend depends on joint_bands and mu_true

  if(is.null(postY)){
    dcip = t(apply(mcmc_output$mu, 2, quantile, c(0.05/2, 1 - 0.05/2)))

    if(include_joint_bands){
      dcib = credBands(mcmc_output$mu)

      plot(t01, y, type='n', ylim=range(dcib, y, na.rm=TRUE), xlab = 't',
           ylab=expression(paste("y"[t])),
           main = 'Credible Intervals for Conditional Expectation',
           cex.lab = 1.5, cex.main = 2, cex.axis = 1)
      polygon(c(t01, rev(t01)), c(dcib[,2], rev(dcib[,1])), col='gray50', border=NA)
      polygon(c(t01, rev(t01)), c(dcip[,2], rev(dcip[,1])), col='grey', border=NA)


      if(!is.null(mu_true)){

        lines(t01, mu_true, lwd=8, col='black', lty=6)
        lines(t01, y, type='p')
        lines(t01, mu, lwd=8, col = 'cyan')
        legend("topleft",
               legend = c("Joint CI", "Pointwise CI", "Pointwise Mean", "True Mean"),
               pch    = c(15, 15, NA, NA),      # squares for the first two
               pt.cex = 2,
               lty    = c(NA, NA, 1, 6),         # no line for the boxes; lines for the last two
               lwd    = c(NA, NA, 8, 8),
               col    = c("gray50", "grey", "cyan", "black"),
               bty    = "n")

      }else{
        lines(t01, y, type='p')
        lines(t01, mu, lwd=8, col = 'cyan')
        legend("topleft",
               legend = c("Joint CI", "Pointwise CI", "Pointwise Mean"),
               pch    = c(15, 15, NA),      # squares for the first two
               pt.cex = 2,
               lty    = c(NA, NA, 1),         # no line for the boxes; lines for the last two
               lwd    = c(NA, NA, 8),
               col    = c("gray50", "grey", "cyan"),
               bty    = "n")
      }

    }else{
      # Plot
      plot(t01, y, type='n', ylim=range(dcip, y, na.rm=TRUE), xlab = 't',
           ylab=expression(paste("y"[t])), main = 'Credible Intervals for Conditional Expectation',
           cex.lab = 1.5, cex.main = 2, cex.axis = 1)
      polygon(c(t01, rev(t01)), c(dcip[,2], rev(dcip[,1])), col='grey', border=NA)

      if(!is.null(mu_true)){

        lines(t01, mu_true, lwd=8, col='black', lty=6)
        lines(t01, y, type='p')
        lines(t01, mu, lwd=8, col = 'cyan')
        legend("topleft",
               legend = c("Pointwise CI", "Pointwise Mean", "True Mean"),
               pch    = c(15, NA, NA),      # squares for the first two
               pt.cex = 2,
               lty    = c(NA, 1, 6),         # no line for the boxes; lines for the last two
               lwd    = c(NA, 8, 8),
               col    = c("grey", "cyan", "black"),
               bty    = "n")

      }else{

        lines(t01, y, type='p')
        lines(t01, mu, lwd=8, col = 'cyan')
        legend("topleft",
               legend = c("Pointwise CI", "Pointwise Mean"),
               pch    = c(15, NA),      # squares for the first two
               pt.cex = 2,
               lty    = c(NA, 1),         # no line for the boxes; lines for the last two
               lwd    = c(NA, 8),
               col    = c("grey", "cyan"),
               bty    = "n")
      }

    }

  }else{
    dcip = t(apply(postY, 2, quantile, c(0.05/2, 1 - 0.05/2)))

    if(include_joint_bands){

      dcib = credBands(postY)

      plot(t01, y, type='n', ylim=range(dcib, y, na.rm=TRUE), xlab = 't',
           ylab=expression(paste("y"[t])), main = 'Credible Intervals for Posterior Predictions',
           cex.lab = 1.5, cex.main = 2, cex.axis = 1)
      polygon(c(t01, rev(t01)), c(dcib[,2], rev(dcib[,1])), col='gray50', border=NA)
      polygon(c(t01, rev(t01)), c(dcip[,2], rev(dcip[,1])), col='grey', border=NA)


      if(!is.null(mu_true)){

        lines(t01, mu_true, lwd=8, col='black', lty=6)
        lines(t01, y, type='p')
        lines(t01, mu, lwd=8, col = 'cyan')
        legend("topleft",
               legend = c("Joint CI", "Pointwise CI", "Pointwise Mean", "True Mean"),
               pch    = c(15, 15, NA, NA),      # squares for the first two
               pt.cex = 2,
               lty    = c(NA, NA, 1, 6),         # no line for the boxes; lines for the last two
               lwd    = c(NA, NA, 8, 8),
               col    = c("gray50", "grey", "cyan", "black"),
               bty    = "n")

      }else{
        lines(t01, y, type='p')
        lines(t01, mu, lwd=8, col = 'cyan')
        legend("topleft",
               legend = c("Joint CI", "Pointwise CI", "Pointwise Mean"),
               pch    = c(15, 15, NA),      # squares for the first two
               pt.cex = 2,
               lty    = c(NA, NA, 1),         # no line for the boxes; lines for the last two
               lwd    = c(NA, NA, 8),
               col    = c("gray50", "grey", "cyan"),
               bty    = "n")
      }

    }else{
      # Plot
      plot(t01, y, type='n', ylim=range(dcip, y, na.rm=TRUE), xlab = 't',
           ylab=expression(paste("y"[t])), main = 'Credible Intervals for Posterior Predictions',
           cex.lab = 1.5, cex.main = 2, cex.axis = 1)
      polygon(c(t01, rev(t01)), c(dcip[,2], rev(dcip[,1])), col='grey', border=NA)

      if(!is.null(mu_true)){

        lines(t01, mu_true, lwd=8, col='black', lty=6)
        lines(t01, y, type='p')
        lines(t01, mu, lwd=8, col = 'cyan')
        legend("topleft",
               legend = c("Pointwise CI", "Pointwise Mean", "True Mean"),
               pch    = c(15, NA, NA),      # squares for the first two
               pt.cex = 2,
               lty    = c(NA, 1, 6),         # no line for the boxes; lines for the last two
               lwd    = c(NA, 8, 8),
               col    = c("grey", "cyan", "black"),
               bty    = "n")

      }else{

        lines(t01, y, type='p')
        lines(t01, mu, lwd=8, col = 'cyan')
        legend("topleft",
               legend = c("Pointwise CI", "Pointwise Mean"),
               pch    = c(15, NA),      # squares for the first two
               pt.cex = 2,
               lty    = c(NA, 1),         # no line for the boxes; lines for the last two
               lwd    = c(NA, 8),
               col    = c("grey", "cyan"),
               bty    = "n")
      }

    }
  }

}

#----------------------------------------------------------------------------
#' Plot the series with change points
#'
#' Plot the time series with distinct segments identified by color.
#'
#' @param mu the \code{T x 1} vector of time series observations (or fitted values)
#' @param cp_inds the \code{n_cp x 1} vector of indices at which a changepoint is identified
#' @export
plot_cp = function(mu, cp_inds){

  dev.new(); par(mfrow=c(1,1), mai = c(1,1,1,1))

  # If no CP's, just plot mu:
  if(length(cp_inds) == 0) return(plot(mu, lwd=8, col =1, type='o'))

  # Assume the CP starts at 1
  if(cp_inds[1]!=1) cp_inds = c(1, cp_inds)

  # Number of changepoints:
  n_cp = length(cp_inds)

  plot(mu, type='n')

  for(j in 1:n_cp) {
    # Indices of CP:
    if(j < n_cp){
      j_ind = cp_inds[j]:(cp_inds[j+1] - 1)
    } else j_ind = cp_inds[length(cp_inds)]:length(mu)

    # Plot in the same color:
    lines(j_ind, mu[j_ind], lwd=8, col =j, type='o')
  }
}
