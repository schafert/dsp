#----------------------------------------------------------------------------
#' Plot the Bayesian trend filtering fitted values
#'
#' Plot the BTF posterior mean of the conditional expectation with posterior credible intervals (pointwise and joint),
#' the observed data, and true curves (if known)
#'
#' @param x an object of class 'dsp' from [dsp_fit()]
#' @param type parameter name; must be included in x$mcmc_output
#' @param true_values (defaults to NULL) the \code{T x 1} vector of the true parameter
#' @param t01 the observation points; if NULL, assume \code{T} equally spaced points from 0 to 1
#' @param include_joint_bands logical; if TRUE, compute simultaneous credible bands (only for `zeta`,`omega`,`yhat`,`mu`)
#' @param ... currently not being used
#'
#' @details
#'
#' The plotting behavior depends on the dimension of the posterior samples stored in `x$mcmc_output[[type]]`:
#' - **1D (scalar parameter):** A density plot is generated using a histogram with overlaid kernel density estimate. The posterior mean and 95% credible interval are annotated, along with the true value if provided.
#' - **2D (vector-valued parameter over time):** A time-series plot is created, showing the posterior mean and 95% pointwise credible intervals. If `include_joint_bands = TRUE` and the parameter is among \code{"omega"}, \code{"mu"}, \code{"yhat"}, or \code{"zeta"}, simultaneous credible bands are also drawn. Optionally, ground truth values (if supplied via `true_values`) are overlaid as orange dots.
#' - **3D (parameter array):** A sequence of time-series plots is drawn, one for each slice of the third dimension (e.g., different components of a multivariate function). Posterior mean, pointwise intervals, joint bands (when applicable), and ground truth are visualized in the same style as the 2D case. The function pauses after each plot, allowing the user to interactively inspect each one.
#' The x-axis values are given by `t01`. If not provided, they default to evenly spaced points in \eqn{[0, 1]}. For parameters with temporal differencing (e.g., `"evol_sigma_t2"`), initial time points used for prior initialization are automatically excluded.
#' If the model includes change point detection (`model = "changepoint"`), and both `omega` and `r` are present in the MCMC output, vertical lines are drawn at the estimated change point locations for plots of `"mu"`, `"yhat"`, or `"omega"`.
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
#' mcmc_output = dsp_fit(y, model_spec = model_spec, nsave = 500, nburn = 500)
#' # Estimated posterior mean vs ground truth
#' plot(mcmc_output, type = "mu", true_values = signal)
#' # Estimated innovation variance vs ground truth for illustration only
#' plot(mcmc_output, type = "obs_sigma_t2", true_values = noise^2)
#'
#' @import coda
#' @importFrom graphics legend plot.new points abline mtext hist legend
#' @importFrom grDevices dev.off devAskNewPage
#' @importFrom stats density predict
#'
#' @method plot dsp
#' @export
plot.dsp = function(x, type ,true_values = NULL, t01 = NULL, include_joint_bands = FALSE, ... ){
  # Time series:
  samples = x$mcmc_output[[type]]
  if(is.null(samples)){
    stop("Must be one of the parameters in x$mcmc_output")
  }
  dimension = dim(samples)
  if(length(dimension) == 2){
    # two dimension -> scatter plot
    # initializing Legend items
    legend_labels <- c("Posterior Mean", "Pointwise CI")
    legend_pch    <- c(NA,15)
    legend_lty    <- c(1,NA)
    legend_lwd    <- c(2,NA)
    legend_col    <- c("steelblue","gray60")

    nT = dimension[2]
    if(is.null(t01)) t01 = seq(0, 1, length.out=nT)

    # Preprocessing for evol_sigma_t2 parameters
    # Dropping first "D" since first D variances are for intializaiton
    if(type %in% c("evol_sigma_t2","zeta_sigma_t2")){
      samples = samples[,-c(1:(x$model_spec$arguments$D))]
      t01 = t01[-c(1:(x$model_spec$arguments$D))]
    }

    # posterior mean
    posterior_mean = colMeans(samples)
    # CI Band
    dcib = credBands(samples)
    dcip = t(apply(samples, 2, quantile, c(0.05/2, 1 - 0.05/2)))

    # setting y limit based on CIs
    ylim = if (include_joint_bands & type %in% c("zeta","omega","yhat","mu"))  range(posterior_mean,dcib,dcip, na.rm = TRUE) else range(posterior_mean,dcip, na.rm = TRUE)

    par(mar = c(4, 4, 1, 9))
    plot.new()
    plot(t01, posterior_mean, type='n', ylim = ylim,
         xlab = '\n',
         ylab='\n')
    if(include_joint_bands & type %in% c("omega","yhat","mu","zeta")){
      polygon(c(t01, rev(t01)), c(dcib[,2], rev(dcib[,1])), col='gray80', border=NA)
      legend_labels <- c(legend_labels, "Joint CI")
      legend_pch    <- c(legend_pch, 15)
      legend_lty    <- c(legend_lty, NA)
      legend_lwd    <- c(legend_lwd, NA)
      legend_col    <- c(legend_col, "gray80")
    }

    polygon(c(t01, rev(t01)), c(dcip[,2], rev(dcip[,1])), col="gray60", border=NA)

    if(!is.null(true_values)){
      if(length(true_values) != length(posterior_mean)){
        dev.off()
        stop("Dimension of true values and the parameters do not match")
      }
      points(t01, true_values, cex = 1.25, col='darkorange',pch = 20)
      legend_labels <- c(legend_labels, "Ground Truth")
      legend_pch    <- c(legend_pch, 20)
      legend_lty    <- c(legend_lty, 1)
      legend_lwd    <- c(legend_lwd, NA)
      legend_col    <- c(legend_col, "darkorange")
    }
    lines(t01, posterior_mean, lwd=5, col = "steelblue")

    # drawing change point?
    draw_cp = x$model_spec$model == "changepoint" &
      all(c("omega","r") %in% names(x$mcmc_output)) &
      (type == "mu" | type == "yhat" | type == "omega")
    if(draw_cp){
      cp_loc = predict(x)
      cp_loc = t01[cp_loc-x$model_spec$arguments$D]
      abline(v = cp_loc,lty = 2,lwd = 2,col = "firebrick")
      legend_labels <- c(legend_labels, "Estimated CP")
      legend_pch    <- c(legend_pch, NA)
      legend_lty    <- c(legend_lty, 2)
      legend_lwd    <- c(legend_lwd, 2)
      legend_col    <- c(legend_col, "firebrick")
    }
    mtext(type, side = 2, line = 2.5,
          cex = 1.25)
    mtext("t",side = 1, line = 2.5,cex = 1.25)

    par(xpd = TRUE)  # turn on drawing outside plot area

    usr <- par("usr")  # get plot coordinate bounds: c(x1, x2, y1, y2)
    x_legend <- usr[2] # a bit to the right
    y_legend <- usr[3]*0.4 + usr[4]*0.6    # vertical center
    graphics::legend(x_legend,y_legend,
                     legend = legend_labels,
                     pch    = legend_pch,
                     lty    = legend_lty,
                     lwd    = legend_lwd,
                     col    = legend_col,
                     pt.cex = 2,
                     bty = "n")
    par(xpd = FALSE)
  }
  else if(is.null(dimension)) {
    # default legend setting
    legend_label = c("Density",
                     "Posterior Mean",
                     "Pointwise CI")
    legend_lty    = c(1,1,2)
    legend_lwd    = c(4,4,2)
    legend_col    = c("black","steelblue","steelblue")

    # one dimension -> density plot
    h = hist(samples,plot = FALSE)
    posterior_mean = mean(samples)
    par(mar = c(4, 4, 1, 9))
    plot.new()
    d <- density(samples)
    ylim <- range(0,h$density,d$y)
    hist(samples, probability = TRUE, col = "gray80", border = "gray80",
         main = "\n", xlab = "\n",ylab = "\n",
         xlim = range(d$x),
         ylim = ylim)
    mtext(type,side = 1, line = 2.5,cex = 1.25)
    mtext("Density",side = 2, line = 2.5,cex = 1.25)
    lines(density(samples), col = "black", lwd = 4)
    abline(v = quantile(samples, c(0.025, 0.975)), col = "steelblue",
           lty = 2,lwd = 2)  # 95% CI
    abline(v = posterior_mean, col = "steelblue", lwd = 4)

    if(!is.null(true_values)) {
      if(length(true_values) != length(posterior_mean)){
        dev.off()
        stop("Dimensions of true values and parameters do not match")
      }
      abline(v = true_values[1],col='darkorange', lwd = 4)
      legend_label = c(legend_label,"Ground Truth")
      legend_lty = c(legend_lty,1)
      legend_lwd = c(legend_lwd,4)
      legend_col = c(legend_col,"darkorange")
    }
    par(xpd = TRUE)
    usr <- par("usr")  # get plot coordinate bounds: c(x1, x2, y1, y2)
    x_legend <- usr[2]
    y_legend <- usr[3]*0.4 + usr[4]*0.6    # vertical center
    graphics::legend(x_legend,
                     y_legend,
                     legend = legend_label,
                     lty    = legend_lty,
                     lwd    = legend_lwd,
                     col    = legend_col,
                     bty = "n")
    par(xpd = FALSE)
  }
  else if(length(dimension) == 3){
    nplots = dimension[3]
    nT = dimension[2]
    if(is.null(t01)) t01 = seq(0, 1, length.out=nT)
    if(type %in% c("evol_sigma_t2","zeta_sigma_t2")){
      samples = samples[,-c(1:(x$model_spec$arguments$D)),]
      t01 = t01[-c(1:(x$model_spec$arguments$D))]
    }
    # setting legend
    legend_labels <- c("Posterior Mean", "Pointwise CI")
    legend_pch    <- c(NA,15)
    legend_lty    <- c(1,NA)
    legend_lwd    <- c(2,NA)
    legend_col    <- c("steelblue","gray60")
    if(include_joint_bands & type %in% c("omega","yhat","mu","beta","zeta")){
      legend_labels <- c(legend_labels, "Joint CI")
      legend_pch    <- c(legend_pch, 15)
      legend_lty    <- c(legend_lty, NA)
      legend_lwd    <- c(legend_lwd, NA)
      legend_col    <- c(legend_col, "gray80")
    }
    if(!is.null(true_values)){
      if(!all(dim(true_values) == c(nT,nplots))){
        stop("Dimension of true values and the parameters do not match.")
      }
      legend_labels <- c(legend_labels, "Ground Truth")
      legend_pch    <- c(legend_pch, 20)
      legend_lty    <- c(legend_lty, 1)
      legend_lwd    <- c(legend_lwd, NA)
      legend_col    <- c(legend_col, "darkorange")
    }
    par(mar = c(4,4,1,9))
    for(i in 1:nplots){
      # posterior mean
      posterior_mean = colMeans(samples[,,i])
      # CI Band
      dcib = credBands(samples[,,i])
      dcip = t(apply(samples[,,i], 2, quantile, c(0.05/2, 1 - 0.05/2)))

      # setting y limit based on CIs
      ylim = if (include_joint_bands & type %in% c("zeta","omega","yhat","mu","beta")) range(posterior_mean,dcib,dcip, na.rm = TRUE) else range(posterior_mean,dcip, na.rm = TRUE)
      plot(t01, posterior_mean, type='n', ylim = ylim,
           xlab = '\n',
           ylab='\n')
      if(include_joint_bands & type %in% c("omega","yhat","mu","zeta","beta")){
        polygon(c(t01, rev(t01)), c(dcib[,2], rev(dcib[,1])), col='gray80', border=NA)
      }
      polygon(c(t01, rev(t01)), c(dcip[,2], rev(dcip[,1])), col="gray60", border=NA)

      if(!is.null(true_values)){
        points(t01, true_values[,i], cex = 1.25, col='darkorange',pch = 20)
      }
      lines(t01, posterior_mean, lwd=5, col = "steelblue")
      mtext(paste0(type,"_",i), side = 2, line = 2.5,
            cex  = 1.25)
      mtext("t",side = 1, line = 2.5,cex = 1.25)


      par(xpd = TRUE)  # turn on drawing outside plot area

      usr <- par("usr")  # get plot coordinate bounds: c(x1, x2, y1, y2)
      x_legend <- usr[2] # a bit to the right
      y_legend <- usr[3]*0.4 + usr[4]*0.6    # vertical center
      graphics::legend(x_legend,y_legend,
                       legend = legend_labels,
                       pch    = legend_pch,
                       lty    = legend_lty,
                       lwd    = legend_lwd,
                       col    = legend_col,
                       pt.cex = 2,
                       bty = "n")
      par(xpd = FALSE)
      devAskNewPage(ask = TRUE)
    }
    #par(ask = FALSE)
    devAskNewPage(ask = FALSE)
    #layout(1)
  }
  else{
    stop("Currently cannot visualize this parameter")
  }
}
