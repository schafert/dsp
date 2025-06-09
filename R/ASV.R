#' MCMC Sampler for Adaptive Stchoastic Volatility (ASV) model
#'
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.
#'
#' @param y the \code{T x 1} vector of time series observations.
#' @param beta the mean of the observed process y. If not provided, they are assumed to be 0.
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1, or D = 2)
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "h" (Log variance)
#' \item "logy2hat" (posterior predictive distribution of log(y^2))
#' \item "sigma2" (Variance, i.e. exp(h))
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#' @param sigma_e
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @export
fit_ASV = function(y,beta = 0,evol_error = "DHS",D = 1,
                   nsave = 1000, nburn = 1000, nskip = 4,
                   mcmc_params = list("h", "logy2hat","sigma2","evol_sigma_t2",
                                      "dhs_phi","dhs_mean"),
                   nugget = FALSE,
                   computeDIC = TRUE,
                   verbose = TRUE){
  nT = length(y)
  y = y - beta
  sample_sd = sd(y)
  # initializing parameters
  if(nugget){
    sParams = init_paramsASV_n(y/sample_sd,evol_error,D)
  }else{
    sParams = init_paramsASV(y/sample_sd,evol_error,D)
  }


  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('h', mcmc_params)) || computeDIC) post_s_mu = array(NA, c(nsave, nT))
  if(!is.na(match('logy2hat', mcmc_params))) post_s_logy2hat = array(NA, c(nsave, nT))
  if(!is.na(match('sigma2', mcmc_params)) || computeDIC) post_sigma2 = array(NA, c(nsave, nT))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_s_evol_sigma_t2 = array(NA, c(nsave, nT))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  if(nugget) post_h_smooth = array(NA,c(nsave,nT))
  if(nugget) sigma2_nugget = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nstot,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
  }
  for(nsi in 1:nstot){
    if((nsi ==10) & verbose){
      pb$tick(10)
    }
    else if( ((nsi%%1000) == 0) & verbose){
      pb$tick(1000)
    }

    if(nugget){
      sParams = fit_paramsASV_n(y/sample_sd,sParams,evol_error,D)
    }else{
      sParams = fit_paramsASV(y/sample_sd,sParams,evol_error,D)
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

        h = sParams$s_mu + log(sample_sd^2)
        h_smooth = sParams$s_mu_sm + log(sample_sd^2)
        sigma2 = exp(h)
        if(!is.na(match('h', mcmc_params)) || computeDIC) post_s_mu[isave,] = h
        if(!is.na(match('logy2hat', mcmc_params))) post_s_logy2hat[isave,] = generate_ly2hat(h,sParams$s_p_error_term)*sample_sd
        if(!is.na(match('sigma2', mcmc_params)) || computeDIC) post_sigma2[isave,] = sigma2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_s_evol_sigma_t2[isave,] = c(sParams$s_evolParams0$sigma_w0^2,
                                                                                         sParams$s_evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = sParams$s_evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = sParams$s_evolParams$dhs_mean
        if(nugget) post_h_smooth[isave,] = h_smooth
        if(nugget) sigma2_nugget[isave] = sParams$s_upevolParams$sigma_wt[1]^2
        post_loglike[isave] = sum(dnorm(y, mean = beta, sd = exp(h/2), log = TRUE))


        # And reset the skip counter:
        skipcount = 0
      }
    }
  }
  if(!is.na(match('h', mcmc_params))) mcmc_output$h = post_s_mu
  if(!is.na(match('logy2hat', mcmc_params))) mcmc_output$logy2hat = post_s_logy2hat
  if(!is.na(match('sigma2', mcmc_params))) mcmc_output$sigma2 = post_sigma2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_s_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean
  if(nugget) mcmc_output$h_smooth = post_h_smooth
  if(nugget) mcmc_output$sigma2_nugget = sigma2_nugget
  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike
  if(computeDIC){
    list_dic_p_d = computeDIC_ASV(y, beta, post_sigma2, post_loglike)
    mcmc_output$DIC = list_dic_p_d$DIC
    mcmc_output$p_d = list_dic_p_d$p_d
  }
  return(mcmc_output);
}
#' Posterior predictive sampler on the transformed y (log(y^2))
#'
#' @param h the log varaince term h
#' @param p_error_term 2 dimensional data frame containing mean and the variance from the 10 componenet
#' Gaussian mixture in Omori et al 2007 paper.
#'
#' @return a vector containing posterior predictive on log(y^2)
#'
generate_ly2hat <- function(h,p_error_term){
  return(h + matrix(rnorm(length(h),
                          mean = p_error_term$mean,
                          sd = sqrt(p_error_term$var)))
  )
}
#' Helper function for initializing parameters for ASV model
#'
#' @param data the \code{T x 1} vector of time series observations.
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1, or D = 2)
#'
#' @return a list containing 4 sets of parameters
#' \itemize{
#' \item s_p_error_term: matrix containing mean and the variance from 10-componenet gaussian mixture (Omori et al. 2007)
#' \item s_mu: a vector containing the posterior sample of log variance h,
#' \item s_evolParams0: a list containing posterior samples of parameters associated with the variance of first D observation of the log variance term, h.
#' \item s_evolParams: a list containing posterior samples parameters associated with the variance of D to the last observations of the log variance temr , h.
#' }
init_paramsASV <- function(data,evol_error,D){
  yoffset = any(data^2 < 10^-16,na.rm = TRUE)*mad(data,na.rm = TRUE)/10^10
  data = log(data^2 + yoffset)
  nT = length(data);
  t01 = seq(0, 1, length.out=nT);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(data)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  data = approxfun(t01, data, rule = 2)(t01)
  loc_obs = t_create_loc(nT, D)
  loc_error = t_create_loc(nT-D,1)

  s_p_error_term = sample_j_wrap(nT,NULL)
  s_mu = sampleBTF(data- s_p_error_term$mean,
                   obs_sigma_t2 = s_p_error_term$var,
                   evol_sigma_t2 = 0.01*rep(1,nT),
                   D = D,
                   loc_obs)
  s_omega = diff(s_mu,differences = D)
  s_mu0 = as.matrix(s_mu[1:D,])
  s_evolParams0 = initEvol0(s_mu0)
  s_evolParams = initEvolParams(s_omega,evol_error)
  return(list(
    s_p_error_term = s_p_error_term,
    s_mu = s_mu,
    s_evolParams0 = s_evolParams0,
    s_evolParams = s_evolParams,
    loc_obs = loc_obs,
    loc_error = loc_error
  ))
}
#' Helper function for Sampling parameters for ASV model
#'
#' @param data the \code{T x 1} vector of time series observations.
#' @param sParams list from the previous run of fit_paramsASV function or init_paramsASV function
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1, or D = 2)
#' @param sigma_e
#'
#' @return a list containing 4 sets of parameters
#' \itemize{
#' \item s_p_error_term: matrix containing mean and the variance from 10-componenet gaussian mixture (Omori et al. 2007)
#' \item s_mu: a vector containing the posterior sample of log variance h,
#' \item s_evolParams0: a list containing posterior samples of parameters associated with the variance of first D observation of the log variance term, h.
#' \item s_evolParams: a list containing posterior samples parameters associated with the variance of D to the last observations of the log variance temr , h.
#' }
fit_paramsASV <- function(data,sParams,evol_error,D){
  yoffset = any(data^2 < 10^-16)*mad(data)/10^10
  data = log(data^2 + yoffset)
  #sigma_e = pi/sqrt(2)
  sigma_e = sd(diff(data,differences = D))
  #sigma_e = 1
  #sigma_e = max(1,sqrt(var(diff(data,differences = D)) - pi^2))
  nT = length(data);
  t01 = seq(0, 1, length.out=nT);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(data)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  data = approxfun(t01, data, rule = 2)(t01)
  #s_p_error_term = sample_jfast(nT,data-sParams$s_mu)
  s_p_error_term = sample_j_wrap(nT,data-sParams$s_mu)
  s_mu = sampleBTF(
    data - s_p_error_term$mean,
    obs_sigma_t2 = s_p_error_term$var,
    evol_sigma_t2 = c(sParams$s_evolParams0$sigma_w0^2,
                      sParams$s_evolParams$sigma_wt^2),
    D = D,
    loc_obs = sParams$loc_obs)
  s_omega = diff(s_mu, differences = D)
  s_mu0 = as.matrix(s_mu[1:D,])
  s_evolParams0 = sampleEvol0(s_mu0, sParams$s_evolParams0)
  s_evolParams = sampleEvolParams(omega = s_omega,
                                  evolParams = sParams$s_evolParams,
                                  sigma_e = sigma_e,
                                  evol_error = evol_error,
                                  loc = sParams$loc_error)
  sParams$s_p_error_term = s_p_error_term
  sParams$s_mu = s_mu
  sParams$s_evolParams0 = s_evolParams0
  sParams$s_evolParams = s_evolParams
  return(sParams)
  #
  # list(s_p_error_term = s_p_error_term,
  #      s_mu = s_mu,
  #      s_evolParams0 = s_evolParams0,
  #      s_evolParams = s_evolParams)
}
#' Helper function for initializing parameters for ASV model  with a nugget effect
#'
#' @param data the \code{T x 1} vector of time series observations.
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1, or D = 2)
#'
#' @return a list containing 4 sets of parameters
#' \itemize{
#' \item s_p_error_term: matrix containing mean and the variance from 10-componenet gaussian mixture (Omori et al. 2007)
#' \item s_mu: a vector containing the posterior sample of log variance h,
#' \item s_evolParams0: a list containing posterior samples of parameters associated with the variance of first D observation of the log variance term, h.
#' \item s_evolParams: a list containing posterior samples parameters associated with the variance of D to the last observations of the log variance temr , h.
#' }
init_paramsASV_n <- function(data,evol_error,D){
  yoffset = any(data^2 < 10^-16,na.rm = TRUE)*mad(data,na.rm = TRUE)/10^10
  data = log(data^2 + yoffset)
  nT = length(data);
  t01 = seq(0, 1, length.out=nT);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(data)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  data = approxfun(t01, data, rule = 2)(t01)
  loc_obs = t_create_loc(nT, D)
  loc_error = t_create_loc(nT-D,1)

  s_p_error_term = sample_j_wrap(nT,NULL)
  # h*
  s_mu = sampleBTF(data- s_p_error_term$mean,
                   obs_sigma_t2 = s_p_error_term$var,
                   evol_sigma_t2 = 0.01*rep(1,nT),
                   D = 0)
  s_upevolParams = initEvolParams(s_mu,evol_error = "NIG")
  # h
  s_mu_sm = sampleBTF(s_mu,
                      obs_sigma_t2 = s_upevolParams$sigma_wt^2,
                      evol_sigma_t2 = 0.01*rep(1,nT),
                      D = D,
                      loc_obs)

  s_omega = diff(s_mu_sm,differences = D)
  s_mu_sm0 = as.matrix(s_mu_sm[1:D,])
  s_evolParams0 = initEvol0(s_mu_sm0)
  s_evolParams = initEvolParams(s_omega,evol_error)


  return(list(
    s_p_error_term = s_p_error_term,
    s_mu = s_mu,
    s_mu_sm = s_mu_sm,
    s_evolParams0 = s_evolParams0,
    s_evolParams = s_evolParams,
    s_upevolParams = s_upevolParams,
    loc_obs = loc_obs,
    loc_error = loc_error
  ))
}
#' Helper function for Sampling parameters for ASV model with a nugget Effect
#'
#' @param data the \code{T x 1} vector of time series observations.
#' @param sParams list from the previous run of fit_paramsASV function or init_paramsASV function
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1, or D = 2)
#'
#' @return a list containing 4 sets of parameters
#' \itemize{
#' \item s_p_error_term: matrix containing mean and the variance from 10-componenet gaussian mixture (Omori et al. 2007)
#' \item s_mu: a vector containing the posterior sample of log variance h,
#' \item s_evolParams0: a list containing posterior samples of parameters associated with the variance of first D observation of the log variance term, h.
#' \item s_evolParams: a list containing posterior samples parameters associated with the variance of D to the last observations of the log variance temr , h.
#' }
fit_paramsASV_n <- function(data,sParams,evol_error,D){
  yoffset = any(data^2 < 10^-16)*mad(data)/10^10
  data = log(data^2 + yoffset)
  #sigma_e = pi/sqrt(2)
  #sigma_e = sd(diff(data,differences = D))
  sigma_e = 1
  #sigma_e = max(1,sqrt(var(diff(data,differences = D)) - pi^2))
  nT = length(data);
  t01 = seq(0, 1, length.out=nT);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(data)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  data = approxfun(t01, data, rule = 2)(t01)
  #s_p_error_term = sample_jfast(nT,data-sParams$s_mu)

  ####################################
  s_p_error_term = sample_j_wrap(nT,data - sParams$s_mu )
  s_mu = sampleBTF(
    y = data - s_p_error_term$mean,
    obs_sigma_t2 = s_p_error_term$var,
    evol_sigma_t2 = sParams$s_upevolParams$sigma_wt^2,
    prior_mean = sParams$s_mu_sm,
    D = 0)

  s_mu_sm = sampleBTF(s_mu,
                      obs_sigma_t2 = sParams$s_upevolParams$sigma_wt^2,
                      evol_sigma_t2 = c(sParams$s_evolParams0$sigma_w0^2,
                                        sParams$s_evolParams$sigma_wt^2),
                      D = D,
                      sParams$loc_obs)

  s_omega = diff(s_mu_sm,differences = D)
  s_mu_sm0 = as.matrix(s_mu_sm[1:D,])
  s_evolParams0 = sampleEvol0(s_mu_sm0, sParams$s_evolParams0)
  s_evolParams = sampleEvolParams(omega = s_omega,
                                  evolParams = sParams$s_evolParams,
                                  sigma_e = sigma_e,
                                  evol_error = evol_error,
                                  loc = sParams$loc_error)

  s_upevolParams <- sampleEvolParams(s_mu - s_mu_sm,
                                     sParams$s_upevolParams,
                                     sigma_e = 1,evol_error = "NIG")

  sParams$s_p_error_term = s_p_error_term
  sParams$s_mu = s_mu
  sParams$s_mu_sm = s_mu_sm
  sParams$s_evolParams0 = s_evolParams0
  sParams$s_evolParams = s_evolParams
  sParams$s_upevolParams = s_upevolParams
  return(sParams)
}
#' Function for calculating DIC and Pb (Bayesian measures of model complexity and fit by Spiegelhalter et al. 2002)
#'
#' @param y the \code{T x 1} vector of time series observations.
#' @param beta the known mean of the process. 0 by default.
#' @param post_sigma2 posterior samples of the variance, i.e. exp(h)
#' @param post_loglike log likelihood based on the posterior sample.
#'
#' @return a list containing DIC and p_d. Two options for estimating both DIC and p_d, which are both included.
#'
computeDIC_ASV <- function(y,beta,post_sigma2,post_loglike){
  # Bayesian measures of model complexity and fit
  # Log-likelihood evaluated at posterior means:
  loglike_hat = sum(dnorm(y,
                          mean = beta,
                          sd = colMeans(sqrt(post_sigma2)),
                          log = TRUE))

  # Effective number of parameters (Note: two options)
  p_d = c(2*(loglike_hat - mean(post_loglike)),
          2*var(post_loglike))
  # DIC:
  DIC = -2*loglike_hat + 2*p_d
  return(list(DIC = DIC, p_d = p_d))
}

