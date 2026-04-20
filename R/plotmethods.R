#----------------------------------------------------------------------------
#' Plot the Bayesian trend filtering fitted values
#'
#' Plot the BTF posterior mean of the conditional expectation with posterior credible intervals (pointwise and joint),
#' the observed data, and true curves (if known)
#'
#' @param x an object of class `dsp` from [dsp_fit()].
#' @param type character string giving the parameter name to visualize; must be one of the entries in `x$mcmc_output`.
#' @param true_values optional ground-truth values to overlay on the plot. For scalar parameters, this should be a length-1 numeric value; for time-varying parameters, a `T x 1` vector; and for multi-parameter time-varying quantities, a `T x p` matrix matching the plotted parameter dimensions.
#' @param y_obs optional vector of observed data point of length T. Only for `2`-dimensional parameters.
#' @param t01 optional vector of observation points. If `NULL`, the function assumes `T` equally spaced points on `[0,1]`.
#' @param include_joint_bands logical; if `TRUE`, include simultaneous credible bands in addition to pointwise credible intervals when available. Joint bands are currently supported only for time-varying parameters: `zeta`, `omega`, `yhat`, and `mu` for 2D outputs, and `zeta`, `omega`, `yhat`, `mu`, and `beta` for 3D outputs.
#' @param alpha numeric credibility level used to construct posterior intervals and bands. Defaults to `0.05`, corresponding to 95% intervals/bands.
#' @param xlab optional x-axis label. If `NULL`, defaults to `"t"`.
#' @param ylab optional y-axis label. If `NULL`, defaults to `type`.
#' @param main optional plot title. For multi-panel plots, this may be either a single title or a vector of titles of length equal to the number of panels.
#' @param xlim optional x-axis limits passed to the plotting routine.
#' @param ylim optional y-axis limits passed to the plotting routine.
#' @param mar optional numeric vector of length 4 giving plot margins, passed to [graphics::par()].
#' @param par_args optional named list of additional graphical parameters passed to [graphics::par()], such as `cex.axis`, `cex.lab`, `cex.main`, `las`, or `mgp`.
#' @param legend logical; if `TRUE`, add a legend to the plot.
#' @param legend_cex numeric scaling factor for legend text size.
#' @param legend_pt_cex numeric scaling factor for legend symbol size.
#' @param nr optional number of rows in the plotting layout for multi-panel (`3`-dimensional) parameters.
#' @param nc optional number of columns in the plotting layout for multi-panel (`3`-dimensional) parameters.
#' @param cp_thres (default 0.5) cutoff proportion for percentage of posterior samples exceeding the threshold needed to label a changepoint
#' @param ... additional graphical arguments passed to the main plotting call, such as [graphics::plot()] or [graphics::hist()] depending on the parameter dimension.
#'
#' @details
#' The plotting behavior depends on the dimension of the posterior samples stored in `x$mcmc_output[[type]]`:
#'
#' \itemize{
#'   \item \strong{1D (scalar parameter):} A density-style summary is produced using a histogram with an overlaid kernel density estimate. The posterior mean and equal-tailed \eqn{(1-\alpha)100\%} credible interval are marked, and the true value is added if supplied through `true_values`.
#'
#'   \item \strong{2D (time-varying parameter):} A time-series plot is produced showing the posterior mean together with equal-tailed pointwise credible intervals. If `include_joint_bands = TRUE` and the parameter is one of `zeta`, `omega`, `yhat`, or `mu`, simultaneous credible bands are also displayed. If `true_values` is provided, the ground truth is overlaid.
#'
#'   \item \strong{3D (multi-parameter time-varying quantity):} A multi-panel collection of time-series plots is produced, one panel for each slice of the third dimension. Each panel shows the posterior mean, pointwise credible intervals, optional joint bands when supported, and optional ground-truth values. The panel layout is controlled by `nr` and `nc`; if omitted, a near-square layout is chosen automatically.
#' }
#'
#' Axis labels, titles, plotting limits, margins, legend display, and additional graphical settings can be customized through `xlab`, `ylab`, `main`, `xlim`, `ylim`, `mar`, `legend`, `legend_cex`, `legend_pt_cex`, and `par_args`.
#'
#' The x-axis values are taken from `t01`. If `t01` is not supplied, evenly spaced points on \eqn{[0,1]} are used. For differenced variance parameters such as `"evol_sigma_t2"` and `"zeta_sigma_t2"`, the initial time points associated with prior initialization are automatically removed before plotting.
#'
#' For fitted changepoint models, changepoint annotations may be added when supported by the plotted parameter and the corresponding latent components are present in the MCMC output.
#'
#' @returns No return value, called for side effects
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
#' @references
#' Krivobokova, T., Kneib, T., and Claeskens, G. (2010).
#' Simultaneous confidence bands for penalized spline estimators.
#' \emph{Journal of the American Statistical Association}, \strong{105}(490), 852--863.
#' \doi{10.1198/jasa.2010.tm09165}
#'
#' @method plot dsp
#' @export

plot.dsp <- function(
  x, type, true_values = NULL, t01 = NULL, y_obs = NULL,
  include_joint_bands = FALSE, alpha = 0.05,
  xlab = NULL, ylab = NULL, main = NULL,
  xlim = NULL, ylim = NULL,
  mar = NULL,
  par_args = list(),      # e.g. list(cex.axis = 0.9, las = 1, mgp = c(2, 0.7, 0))
  legend = TRUE,
  legend_cex = 1,
  legend_pt_cex = 2,
  nr = NULL,
  nc = NULL,
  cp_thres = 0.5,
  ...
){
  # Time series:
  mean_color = "dodgerblue3"
  samples = x$mcmc_output[[type]]
  if(is.null(samples)){
    stop("Must be one of the parameters in x$mcmc_output")
  }
  dimension = dim(samples)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # Apply par settings in one place (default + user override)
  if (!is.null(mar)) par(mar = mar) else par(mar = c(4, 4, 1, 9))
  if (length(par_args)) do.call(par, par_args)

  plot_args <- list(...)  # optional pass-through for plot/points/lines/hist etc.

  # Helper defaults
  default_xlab <- "t"
  default_ylab <- type
  xlab_use <- if (is.null(xlab)) default_xlab else xlab
  ylab_use <- if (is.null(ylab)) default_ylab else ylab
  main_use <- if (is.null(main)) "" else main

  # example: 2D branch
  if(length(dimension) == 2){
    # two dimension -> scatter plot
    # initializing Legend items
    nT = dimension[2]
    if(is.null(t01)) t01 = seq(0, 1, length.out=nT)

    # Preprocessing for evol_sigma_t2 parameters
    # Dropping first "D" since first D variances are for intializaiton
    if(type %in% c("evol_sigma_t2","zeta_sigma_t2")){
      samples = samples[,-c(1:(x$model_spec$arguments$D))]
      t01 = t01[-c(1:(x$model_spec$arguments$D))]
    }

    posterior_mean <- colMeans(samples)
    dcip <- t(apply(samples, 2, quantile, c(alpha/2, 1-alpha/2)))
    dcib = credBands(samples,alpha)

    ylim_use <- if (is.null(ylim)) {
      if (include_joint_bands && type %in% c("zeta", "omega", "yhat", "mu")) {
        range(posterior_mean, dcib, dcip, na.rm = TRUE)
      } else {
        range(posterior_mean, dcip, na.rm = TRUE)
      }
    } else {
      ylim
    }
    xlim_use <- if (is.null(xlim)) range(t01, na.rm = TRUE) else xlim

    legend_labels <- c("Posterior Mean", "Pointwise CI")
    legend_pch <- c(NA, 15)
    legend_lty <- c(1, NA)
    legend_lwd <- c(3, NA)
    legend_col <- c(mean_color, "gray75")

    # changepoint block
    draw_cp <- x$model_spec$model == "changepoint" && all(c("omega", "r") %in% names(x$mcmc_output)) && type %in% c("mu", "yhat", "omega")
    if(draw_cp){
      legend_labels <- c(legend_labels, "Estimated CP")
      legend_pch <- c(legend_pch, NA)
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 2)
      legend_col <- c(legend_col, "firebrick")
      cp_loc = predict(x,cp_thres = cp_thres)
      cp_loc = t01[cp_loc - x$model_spec$arguments$D]
    }
    # build plot() args cleanly
    base_plot_args <- c(
      list(
        x = t01, y = posterior_mean, type = "n",
        xlab = xlab_use, ylab = ylab_use, main = main_use,
        xlim = xlim_use, ylim = ylim_use
      ),
      plot_args
    )
    do.call(plot, base_plot_args)
    if(include_joint_bands && type %in% c("omega", "yhat", "mu", "zeta")){
      polygon(
        c(t01, rev(t01)),
        c(dcib[, 2], rev(dcib[, 1])),
        col = "gray85", border = NA
      )
      legend_labels <- c(legend_labels, "Joint CI")
      legend_pch <- c(legend_pch, 15)
      legend_lty <- c(legend_lty, NA)
      legend_lwd <- c(legend_lwd, NA)
      legend_col <- c(legend_col, "gray85")
    }

    polygon(
      c(t01, rev(t01)),
      c(dcip[, 2], rev(dcip[, 1])),
      col = "gray75", border = NA
    )
    if(draw_cp){abline(v = cp_loc, lty = 2, lwd = legend_lwd, col = "firebrick")}
    lines(t01, posterior_mean, lwd = 3, col = mean_color)
    if (!is.null(true_values)) {
      points(t01, true_values, pch = 20, col = "darkorange",cex = legend_pt_cex)
      legend_labels <- c(legend_labels, "Ground Truth")
      legend_pch <- c(legend_pch, 20)
      legend_lty <- c(legend_lty, NA)
      legend_lwd <- c(legend_lwd, NA)
      legend_col <- c(legend_col, "darkorange")
    }

    if (!is.null(y_obs)) {
      if (length(y_obs) != length(t01)) {
        stop("Length of y_obs must match the plotted time dimension.")
      }
      points(t01, y_obs, pch = 16, col = "black", cex = 0.7)
    }

    if (legend) {
      usr <- par("usr")
      par(xpd = TRUE)
      on.exit(par(xpd = FALSE), add = TRUE)
      graphics::legend(
        usr[2], usr[3] * 0.4 + usr[4] * 0.6,
        legend = legend_labels,
        pch = legend_pch,
        lty = legend_lty,
        lwd = legend_lwd,
        col = legend_col,
        pt.cex = legend_pt_cex,
        cex = legend_cex,
        bty = "n"
      )
    }
  }else if(is.null(dimension) || length(dimension) == 1) {
    # default legend setting
    legend_labels = c("Density",
                     "Posterior Mean",
                     "Pointwise CI")
    legend_lty    = c(1,1,2)
    legend_lwd    = c(4,4,2)
    legend_col    = c("black",mean_color,mean_color)

    # one dimension -> density plot
    h = hist(samples,plot = FALSE)
    posterior_mean = mean(samples)
    d <- density(samples)

    xlim_use <- if (is.null(xlim)) range(d$x, na.rm = TRUE) else xlim
    ylim_use <- if (is.null(ylim)) range(0, h$density, d$y, na.rm = TRUE) else ylim
    base_hist_args <- c(
      list(
        x = samples,
        probability = TRUE,
        col = "gray85",
        border = "gray85",
        main = main_use,
        xlab = xlab_use,
        ylab = ylab_use,
        xlim = xlim_use,
        ylim = ylim_use
      ),
      plot_args
    )
    do.call(hist, base_hist_args)
    lines(d, col = "black", lwd = 4)
    abline(v = quantile(samples, c(alpha/2, 1 - alpha/2)), col = mean_color, lty = 2, lwd = 2)
    abline(v = posterior_mean, col = mean_color, lwd = 4)

    if (!is.null(true_values)) {
      if (length(true_values) != 1) {
        stop("Dimensions of true values and parameters do not match")
      }
      abline(v = true_values[1], col = "darkorange", lwd = 4)
      legend_labels <- c(legend_labels, "Ground Truth")
      legend_lty <- c(legend_lty, 1)
      legend_lwd <- c(legend_lwd, 4)
      legend_col <- c(legend_col, "darkorange")
    }

    if (legend) {
      usr <- par("usr")
      par(xpd = TRUE)
      on.exit(par(xpd = FALSE), add = TRUE)
      graphics::legend(
        usr[2], usr[3] * 0.4 + usr[4] * 0.6,
        legend = legend_labels,
        lty = legend_lty,
        lwd = legend_lwd,
        col = legend_col,
        pt.cex = legend_pt_cex,
        cex = legend_cex,
        bty = "n"
      )
    }
  }
  else if(length(dimension) == 3){
    nplots = dimension[3]
    nT = dimension[2]
    if(is.null(t01)) t01 = seq(0, 1, length.out=nT)
    # Preprocessing for evol_sigma_t2 parameters
    # Dropping first "D" since first D variances are for intializaiton
    if(type %in% c("evol_sigma_t2", "zeta_sigma_t2")){
      samples = samples[, -c(1:(x$model_spec$arguments$D)), ]
      t01 = t01[-c(1:(x$model_spec$arguments$D))]
    }

    if(!is.null(true_values)){
      if(!all(dim(true_values) == c(length(t01), nplots))){
        stop("Dimension of true values and the parameters do not match.")
      }
    }
    # setting legend
    legend_labels <- c("Posterior Mean", "Pointwise CI")
    legend_pch <- c(NA, 15)
    legend_lty <- c(1, NA)
    legend_lwd <- c(3, NA)
    legend_col <- c(mean_color, "gray75")

    if(include_joint_bands && type %in% c("omega", "yhat", "mu", "beta", "zeta")){
      legend_labels <- c(legend_labels, "Joint CI")
      legend_pch <- c(legend_pch, 15)
      legend_lty <- c(legend_lty, NA)
      legend_lwd <- c(legend_lwd, NA)
      legend_col <- c(legend_col, "gray85")
    }

    if(!is.null(true_values)){
      legend_labels <- c(legend_labels, "Ground Truth")
      legend_pch <- c(legend_pch, 20)
      legend_lty <- c(legend_lty, NA)
      legend_lwd <- c(legend_lwd, NA)
      legend_col <- c(legend_col, "darkorange")
    }

    # layout
    if(is.null(nr) && is.null(nc)){
      nr_use <- ceiling(sqrt(nplots))
      nc_use <- ceiling(nplots / nr_use)
    } else if(is.null(nr)) {
      nc_use <- nc
      nr_use <- ceiling(nplots / nc_use)
    } else if(is.null(nc)) {
      nr_use <- nr
      nc_use <- ceiling(nplots / nr_use)
    } else {
      nr_use <- nr
      nc_use <- nc
    }

    if(nr_use * nc_use < nplots){
      stop("nr * nc must be at least the number of panels.")
    }

    par_args$mfrow <- c(nr_use, nc_use)
    do.call(par, par_args)

    for(i in 1:nplots){
      posterior_mean <- colMeans(samples[, , i])
      dcip <- t(apply(samples[, , i], 2, quantile, c(alpha/2, 1 - alpha/2)))
      dcib <- credBands(samples[, , i], alpha)

      ylim_use <- if (is.null(ylim)) {
        if (include_joint_bands && type %in% c("zeta", "omega", "yhat", "mu", "beta")) {
          range(posterior_mean, dcib, dcip, na.rm = TRUE)
        } else {
          range(posterior_mean, dcip, na.rm = TRUE)
        }
      } else {
        ylim
      }

      xlim_use <- if (is.null(xlim)) range(t01, na.rm = TRUE) else xlim
      main_use_i <- if (is.null(main) || length(main) != nplots) paste0(type, "_", i) else main[i]

      base_plot_args <- c(
        list(
          x = t01, y = posterior_mean, type = "n",
          xlab = xlab_use, ylab = ylab_use, main = main_use_i,
          xlim = xlim_use, ylim = ylim_use
        ),
        plot_args
      )
      do.call(plot, base_plot_args)

      if(include_joint_bands && type %in% c("omega", "yhat", "mu", "zeta", "beta")){
        polygon(
          c(t01, rev(t01)),
          c(dcib[, 2], rev(dcib[, 1])),
          col = "gray85", border = NA
        )
      }

      polygon(
        c(t01, rev(t01)),
        c(dcip[, 2], rev(dcip[, 1])),
        col = "gray75", border = NA
      )

      if(!is.null(true_values)){
        points(t01, true_values[, i], pch = 20, col = "darkorange",cex = legend_pt_cex)
      }

      lines(t01, posterior_mean, lwd = 3, col = mean_color)

      if (legend) {
        usr <- par("usr")
        par(xpd = TRUE)
        on.exit(par(xpd = FALSE), add = TRUE)
        graphics::legend(
          usr[2], usr[3] * 0.4 + usr[4] * 0.6,
          legend = legend_labels,
          pch = legend_pch,
          lty = legend_lty,
          lwd = legend_lwd,
          col = legend_col,
          pt.cex = legend_pt_cex,
          cex = legend_cex,
          bty = "n"
        )
        par(xpd = FALSE)
      }

    }
  }else{
    stop("Currently cannot visualize this parameter")
  }
}
