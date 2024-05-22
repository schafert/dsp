#' MCMC Sampler for Dynamic Shrinkage with Changepoint
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on zeroth (D = 0),
#' first (D = 1), or second (D = 2) differences of the conditional expectation.
#' The user has the option to utilize threshold shrinkage to identify changepoints
#' for D = 1 or 2.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Note that only 'DHS' works for threshold shrinkage with changepoints.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param cp flag indicator to determine whether to use threshold shrinkage with changepoints.
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 0, D = 1, or D = 2)
#' @param useObsSV logical; if TRUE, include a (normal) stochastic volatility model
#' for the observation error variance
#' @param useAnom logical; if TRUE, include an anomaly component in the observation equation
#' (only for threshold shrinkage with changepoints.)
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#' @param cp_thres Percentage of posterior samples needed to declare a changepoint
#' @param return_full_samples logical; if TRUE, return full MCMC samples for desired parameters;
#' if FALSE, return the posterior meean for desired parameters
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' if threshold shrinkage with changepoints is used, also return detected changepoint locations
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#' \dontrun{
#' # Example 1: Change in mean with stochastic volatility
#' signal = c(rep(0, 50), rep(10, 50))
#' noise = rep(1, 100)
#' noise_var = rep(1, 100)
#' for (k in 2:100){
#'   noise_var[k] = exp(0.9*log(noise_var[k-1]) + rnorm(1, 0, 0.5))
#'   noise[k] = rnorm(1, 0, sqrt(noise_var[k])) }
#'
#' y = signal + noise
#' mcmc_output = dsp_cp(y, cp=TRUE, mcmc_params = list('yhat', 'mu', "omega", "r"))
#' cp = mcmc_output$cp
#' print(paste0('Changepoint Locations: ', cp))
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)
#'
#' # Example 2: Change in mean with Outliers
#' signal = c(rep(0, 100), rep(5, 100))
#' y = c(rep(0, 100), rep(5, 100)) + rnorm(200)
#' y[50] = 10
#' y[150] = -10
#'
#' mcmc_output = dsp_cp(y, cp=TRUE, useAnom = TRUE, mcmc_params = list('yhat', 'mu', "omega", "r"))
#' cp = mcmc_output$cp
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)
#'
#' # Example 3: Change in linear trend
#' signal = c(seq(1, 50), seq(51, 2))
#' y = c(seq(1, 50), seq(51, 2)) + rnorm(100)
#'
#' mcmc_output = dsp_cp(y, cp=TRUE, D=2, mcmc_params = list('yhat', 'mu', "omega", "r"))
#' cp = mcmc_output$cp
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = signal)
#' }
#'
#' @export
dsp_cp = function(y, cp = FALSE, evol_error = 'DHS', D = 1, useObsSV = TRUE, useAnom = FALSE,
               nsave = 1000, nburn = 1000, nskip = 4,
               mcmc_params = list("mu", "omega", "r"),
               computeDIC = TRUE,
               verbose = TRUE,
               cp_thres = 0.4,
               return_full_samples = TRUE){
  if(!((evol_error == "DHS") || (evol_error == "HS") || (evol_error == "BL") || (evol_error == "SV") || (evol_error == "NIG"))) stop('Error type must be one of DHS, HS, BL, SV, or NIG')
  if(!((D == 0) || (D == 1) || (D == 2))) stop('D must be 0, 1 or 2')

  if (cp == TRUE) {
    if (evol_error == 'DHS' & (D == 1 || D == 2)) {
      # Add in omega and r for finding changepoints
      if (is.na(match('omega', mcmc_params))) {
        mcmc_params = append(mcmc_params, list('omega'))
      }
      if (is.na(match('r', mcmc_params))) {
        mcmc_params = append(mcmc_params, list('r'))
      }
      mcmc_output = abco(y, D = D, useObsSV = useObsSV, useAnom = useAnom, nsave = nsave,
                         nburn = nburn, nskip = nskip, mcmc_params = mcmc_params, verbose = verbose, cp_thres = cp_thres)
    }
  } else {
    mcmc_output = btf(y, evol_error = evol_error, D = D, useObsSV = useObsSV, nsave = nsave, nburn = nburn, nskip = nskip,
              mcmc_params = mcmc_params, computeDIC = computeDIC, verbose = verbose)
  }

  if (!return_full_samples){
    for (nm in names(mcmc_output)){
      if (!is.na(match(nm, mcmc_params))) {
        mcmc_output[[nm]] = colMeans(as.matrix(mcmc_output[[nm]]))
      }
    }
  }
  return (mcmc_output);
}
#' MCMC Sampler for Bayesian Trend Filtering
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on zeroth (D = 0),
#' first (D = 1), or second (D = 2) differences of the conditional expectation.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior; the default), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 0, D = 1, or D = 2)
#' @param useObsSV logical; if TRUE, include a (normal) stochastic volatility model
#' for the observation error variance
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burnin)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#' \dontrun{
#' # Example 1: Bumps Data
#' simdata = simUnivariate(signalName = "bumps", T = 128, RSNR = 7, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf(y, D=1)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#' # Example 2: Doppler Data; longer series, more noise
#' simdata = simUnivariate(signalName = "doppler", T = 500, RSNR = 5, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf(y)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#' # And examine the AR(1) parameters for the log-volatility w/ traceplots:
#' plot(as.ts(out$dhs_phi)) # AR(1) coefficient
#' plot(as.ts(out$dhs_mean)) # Unconditional mean
#'
#'# Example 3: Blocks data (locally constant)
#' simdata = simUnivariate(signalName = "blocks", T = 1000, RSNR = 3, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf(y, D = 1) # try D = 1 to approximate the locally constant behavior
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#' }
#'
#' @export
btf = function(y, evol_error = 'DHS', D = 2, useObsSV = FALSE,
               nsave = 1000, nburn = 1000, nskip = 4,
               mcmc_params = list("mu", "yhat","evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
               computeDIC = TRUE,
               verbose = TRUE){

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # For D = 0, return special case:
  if(D == 0){
    return(btf0(y = y, evol_error = evol_error, useObsSV = useObsSV,
                nsave = nsave, nburn = nburn, nskip = nskip,
                mcmc_params = mcmc_params,
                computeDIC = computeDIC,
                verbose = verbose))
  }

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Initialize bandsparse locations
  loc = t_create_loc(length(y)-D, 1)
  loc_obs = t_create_loc(length(y), D)

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Compute the Cholesky term (uses random variances for a more conservative sparsity pattern)
  chol0 = NULL

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = D, loc_obs = loc_obs, chol0)

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  # And the initial states:
  mu0 = as.matrix(mu[1:D,])

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(mu0)

  # SV parameters, if necessary:
  if(useObsSV) {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}

  # For HS MCMC comparisons:
  # evolParams$dhs_phi = 0

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the states:
    mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2), D = D, loc_obs, chol0)

    # Compute the evolution errors:
    omega = diff(mu, differences = D)

    # And the initial states:
    mu0 = as.matrix(mu[1:D,])

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(mu0, evolParams0, A = 1)

    # Sample the (observation and evolution) variances and associated parameters:
    if(useObsSV){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, 1/sqrt(T), evol_error, loc)

      # Observation error variance + params:
      svParams = sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt

    } else {
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(T), evol_error, loc)

      if(evol_error == 'DHS') {
        sigma_e = uni.slice(sigma_e, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum(evolParams$xiLambda)))
      if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum((omega/evolParams$tau_j)^2)/2))
      if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

      # Replicate for coding convenience:
      sigma_et = rep(sigma_e, T)
    }

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
#' MCMC Sampler for Bayesian Trend Filtering: D = 0
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on the conditional expectation.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param useObsSV logical; if TRUE, include a (normal) stochastic volatility model
#' for the observation error variance
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burnin)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
btf0 = function(y, evol_error = 'DHS', useObsSV = FALSE,
                nsave = 1000, nburn = 1000, nskip = 4,
                mcmc_params = list("mu", "yhat","evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                computeDIC = TRUE,
                verbose = TRUE){

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);
  loc = t_create_loc(length(y), 0)

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = 0)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega = mu, evol_error = evol_error)

  # SV parameters, if necessary:
  if(useObsSV) {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the states:
    mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = evolParams$sigma_wt^2, D = 0)

    # Sample the (observation and evolution) variances and associated parameters:
    if(useObsSV){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega = mu, evolParams, 1/sqrt(T), evol_error, loc)

      # Observation error variance + params:
      svParams = sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt

    } else {
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega = mu, evolParams, sigma_e/sqrt(T), evol_error, loc)

      # Sample the observation error SD:
      if(evol_error == 'DHS') {
        sigma_e = uni.slice(sigma_e, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum(evolParams$xiLambda)))
      if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum((mu/evolParams$tau_j)^2)/2))
      if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

      # Replicate for coding convenience:
      sigma_et = rep(sigma_e, T)
    }

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] =  evolParams$sigma_wt^2
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
