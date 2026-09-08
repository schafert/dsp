#----------------------------------------------------------------------------
#' Generate univariate signals of different type
#'
#' Using code from the archived \code{wmtsa} package
#'
#' @param name character string of name of the test wavelet signal to be generated; one of "dirac", "kronecker", "heavisine", "bumps", "blocks",
#' "doppler", "ramp", "cusp", "crease", "sing", "hisine",
#' "losine", "linchirp", "twochirp", "quadchirp",
#' "mishmash1", "mishmash2", "mishmash3", "levelshift",
#' "jumpsine", "gauss", "patches",
#' "linear", "quadratic", "cubic";
#' @param n length of the series; defaults to 1024 points; increasing n infills the time series
#' @param snr desired signal-to-noise ratio; default \code{Inf} corresponds to 0 noise
#'
#' @returns A numeric vector the same length as `n`.
#'
#' @examples
#' nms <- c("blocks", "linchirp", "mishmash1", "bumps")
#' z <- lapply(nms, simUnivariate)
#'
#' @export

simUnivariate <- function(name, n=1024, snr=Inf)
{

  .wave.demo.signals <- c("dirac", "kronecker", "heavisine", "bumps", "blocks",
                            "doppler", "ramp", "cusp", "crease", "sing", "hisine",
                            "losine", "linchirp", "twochirp", "quadchirp",
                            "mishmash1", "mishmash2", "mishmash3", "levelshift",
                            "jumpsine", "gauss", "patches",
                            "linear", "quadratic", "cubic")

  x <- (0:(n-1.))/n
  z <- switch(name,
              dirac=n*(x == floor(.37*n)/n),
              kronecker=(x == floor(.37*n)/n),
              heavisine=4*sin(4*pi*x)-sign(x-.3)-sign(.72-x),
              bumps={
                pos <- c(.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81)
                hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
                wth <- c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008,.005)
                y <- rep(0, n)
                for(j in 1:length(pos)) y <- y+hgt[j]/(1+abs((x-pos[j]))/wth[j])^4
                y
              },
              blocks={
                pos <- c(.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81)
                hgt <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1,2.1, -4.2)
                y <- rep(0, n)
                for(j in 1:length(pos)) y <- y+(1+sign(x-pos[j]))*hgt[j]/2
                y
              },
              doppler=sqrt(x*(1-x))*sin((2*pi*1.05)/(x+.05)),
              ramp=x-(x >= .37),
              cusp=sqrt(abs(x-.37)),
              crease=exp(-4*abs(x-.5)),
              sing=1/abs(x-(floor(n*.37)+.5)/n),
              hisine=sin(pi*n*.6902*x),
              midsine=sin(pi*n*.3333*x),
              losine=sin(pi*n*.03*x),
              linchirp=sin(.125*pi*n*x^2),
              twochirp=sin(pi*n*x^2) + sin((pi/3)*n*x^2),
              quadchirp=sin((pi/3)*n*x^3),
              # QuadChirp + LinChirp + HiSine
              mishmash1=sin((pi/3)*n*x^3) + sin(pi*n*.6902*x) + sin(pi*n*.125*x^2),
              # QuadChirp + LinChirp + HiSine + Bumps
              mishmash2={		# wernersorrows
                y   <- sin(pi*(n/2)*x^3)+sin(pi*n*.6902*x)+sin(pi*n*x^2)
                pos <- c(.1, .13, .15, .23, .25, .40, .44, .65, .76, .78, .81)
                hgt <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)
                wth <- c(.005, .005, .006, .01, .01, .03, .01, .01, .005, .008,.005)
                for(j in 1:length(pos)) y <- y + hgt[j]/(1+abs((x-pos[j])/wth[j]))^4
                y
              },
              # QuadChirp + MidSine + LoSine + Sing/200.
              mishmash3=sin((pi/3)*n*x^3) + sin(pi*n*.3333*x) + sin(pi*n*.03*x) +
                (1/abs(x-(floor(n*.37)+.5)/n))/(200.*n/512.),
              gauss=dnorm(x, .3, .025),
              jumpsine=10.*(sin(4*pi*x) + as.numeric(x >= 0.625 & x < 0.875)),
              levelshift=as.numeric(x >= 0.25 & x < 0.39),
              patches={
                if(n<16) stop("n must be >= 16 to generate patches\n")
                J <- logb(n, base=2)
                y <- rep(0., n)
                for(j in 0:(J-4.)) y[(1:2^j)+3.*2.^(j+2.)] <- 1.
                y
              },
              linear=2.*x-1.,
              quadratic=4. * (1. - x) * x,
              cubic=64. * x * (x - 1.) * (x - .5) / 3.,
              stop("Unknown signal name.  Allowable names are:\n",
                   paste(.wave.demo.signals, collapse = ", ")))

  if (snr > 0)
    z <- z + rnorm(n) * sqrt(var(z)) / snr

  z
}

#----------------------------------------------------------------------------
#' Simulate noisy observations from a dynamic regression model
#'
#' Simulates data from a time series regression with dynamic regression coefficients.
#' The dynamic regression coefficients are simulated as a Gaussian random walk,
#' where jumps occur with a pre-specified probability \code{sparsity}.
#' The coefficients are initialized by a N(0,1) simulation.
#'
#' @param nT number of time points
#' @param p number of predictors (total)
#' @param p_0 number of true zero regression terms
#' @param sparsity the probability of a jump
#' (i.e., a change in the dynamic regression coefficient)
#' @param RSNR root-signal-to-noise ratio
#' @param ar1 the AR(1) coefficient for the predictors X; default is zero for iid N(0,1) predictors
#' @param include_plot logical; if TRUE, include a plot of the simulated data and the true curve
#'
#' @return a list containing
#' \itemize{
#' \item the simulated function \code{y}
#' \item the simulated predictors \code{X}
#' \item the simulated dynamic regression coefficients \code{beta_true}
#' \item the true function \code{mu_true}
#' \item the true observation standard deviation \code{sigma_true}
#' }
#'
#'
#' @note The root-signal-to-noise ratio is defined as RSNR = (sd of true function)/(sd of noise).
#' @importFrom stats arima.sim
#' @export

simRegression = function(nT = 200, p = 20, p_0 = 15,
                         sparsity = 0.05, RSNR = 5, ar1 = 0,
                         include_plot = FALSE){

  if(p < p_0) stop('Must have more predictors (p) than true zeros (p_0)')

  # Simulate the predictors: autocorrelated or independent?
  # Either way, use N(0,1) innovations
  if(ar1 == 0){
    X = cbind(1,matrix(rnorm(n = nT*(p-1)), nrow = nT, ncol = p-1))
  } else X = cbind(1,
                   apply(matrix(0, nrow = nT, ncol = p-1), 2, function(x)
                     stats::arima.sim(n = nT, list(ar = ar1), sd = sqrt(1-ar1^2))))

  # Simulate the true regression signals
  beta_true = matrix(0, nrow = nT, ncol = p);

  # Value of intercept:
  beta_true[,1] = 1

  # Now, for the remaining nonzero predictors, simulate as jumps:
  if((p - p_0) > 1){for(j in 2:(p - p_0)){
    # Simulate the paths:
    beta_true[,j] = rnorm(n = 1) +
      cumsum(rnorm(n = nT)*rbinom(n = nT, size = 1, prob = sparsity))
  }}

  # Conditional mean:
  mu_true = rowSums(X*beta_true)

  # Noise SD, based on RSNR (also put in a check for constant/zero functions)
  sigma_true = sd(mu_true)/RSNR; if(sigma_true==0) sigma_true = sqrt(sum(mu_true^2)/nT)/RSNR + 10^-3

  # Observed data:
  y = mu_true + sigma_true*rnorm(nT)

  # Plot?
  if(include_plot) {t = seq(0, 1, length.out=nT); plot(t, y, main = 'Simulated Data and True Curve'); lines(t, mu_true, lwd=8, col='black') }

  # Return the raw data and the true values:
  list(y = y, X = X, beta_true = beta_true, mu_true = mu_true, sigma_true = sigma_true)
}
#----------------------------------------------------------------------------
#' Simulate noisy observations from a dynamic regression model
#'
#' Simulates data from a time series regression with dynamic regression coefficients.
#' The dynamic regression coefficients are selected using the options from the
#' \code{simUnivariate()} function in the \code{wmtsa} package.
#'
#' @param signalNames vector of strings matching the "name" argument in the \code{simUnivariate()} function,
#' e.g. "bumps" or "doppler"
#' @param nT number of points
#' @param RSNR root-signal-to-noise ratio
#' @param p_0 number of true zero regression terms to include
#' @param include_intercept logical; if TRUE, the first column of X is 1's
#' @param scale_all logical; if TRUE, scale all regression coefficients to \[0,1\]
#' @param include_plot logical; if TRUE, include a plot of the simulated data and the true curve
#' @param ar1 the AR(1) coefficient for the predictors X; default is zero for iid N(0,1) predictors
#'
#' @return a list containing
#' \itemize{
#' \item the simulated function \code{y}
#' \item the simulated predictors \code{X}
#' \item the simulated dynamic regression coefficients \code{beta_true}
#' \item the true function \code{mu_true}
#' \item the true observation standard deviation \code{sigma_true}
#' }
#'
#' @note The number of predictors is \code{p = length(signalNames) + p_0}.
#'
#' @note The root-signal-to-noise ratio is defined as RSNR = (sd of true function)/(sd of noise).
#'
#' @keywords internal
simRegression0 = function(signalNames = c("bumps", "blocks"), nT = 200, RSNR = 10, p_0 = 5, include_intercept = TRUE, scale_all = TRUE, include_plot = TRUE, ar1 = 0){

  # True number of signals
  p_true = length(signalNames)

  # Total number of covariates (non-intercept)
  p = p_true + p_0

  # Simulate the true regression signals
  beta_true = matrix(0, nrow = nT, ncol = p)
  for(j in 1:p_true) beta_true[,j] = simUnivariate(signalNames[j], n=nT);
  if(scale_all) beta_true[,1:p_true] = apply(as.matrix(beta_true[,1:p_true]), 2, function(x) (x - min(x))/(max(x) - min(x)))

  # Simulate the predictors: autocorrelated or independent? Either way, use N(0,1) innovations
  if(ar1 == 0){
    X = matrix(rnorm(nT*p), nrow=nT, ncol = p)
  } else X = apply(matrix(0, nrow = nT, ncol = p), 2, function(x) stats::arima.sim(n = nT, list(ar = ar1), sd = sqrt(1-ar1^2)))

  # If we want an intercept, simply replace the first column w/ 1s
  if(include_intercept) X[,1] = matrix(1, nrow = nrow(X), ncol = 1)

  # The true response function:
  mu_true = rowSums(X*beta_true)

  # Noise SD, based on RSNR (also put in a check for constant/zero functions)
  sigma_true = sd(mu_true)/RSNR; if(sigma_true==0) sigma_true = sqrt(sum(mu_true^2)/nT)/RSNR + 10^-3

  # Simulate the data:
  y = mu_true + sigma_true*rnorm(nT)

  # Plot?
  if(include_plot) {t = seq(0, 1, length.out=nT); plot(t, y, main = 'Simulated Data and True Curve'); lines(t, mu_true, lwd=8, col='black') }

  # Return the raw data and the true values:
  list(y = y, X = X, beta_true = beta_true, mu_true = mu_true, sigma_true = sigma_true)
}
