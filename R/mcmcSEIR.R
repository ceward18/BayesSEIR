################################################################################
# General MCMC function
# takes in data and infPeriodSpec and fits specified MCMC algorithm
################################################################################

#' Run the MCMC algorithm to obtain posterior samples for all model parameters.
#'
#' @param dat A \emph{list} providing the necessary data on the epidemic. Must contain 
#' \code{Istar, S0, E0, I0, and N}. If \code{EKnown = TRUE}, must also contain Estar.
#' @param X The design matrix corresponding to the intensity process over epidemic time.
#' @param inits A \emph{list} of initial values for all parameters. Must contain initial
#' values for beta vector with elements for each column of \code{X}. Other parameters
#' needed depend on the model being fit - see Details.
#' @param niter The total number of iterations to be run.
#' @param nburn The number of burn-in iterations to be discarded.
#' @param infPeriodSpec A character indicating which model for the infectious period
#' will be used. 
#' @param priors A \emph{list} of functions to compute the prior probability for each parameter.
#' Functions must take only one argument and must compute the probability on the log scale.
#' @param iddFun If \code{infPeriodSpec = "IDD"} this will specify the functional 
#' form used to describe the IDD transmissibility curve. Default functions are 
#' \code{dgammaIDD}, \code{dlnormIDD}, \code{logitIDD}, and \code{splineIDD}
#' @param maxInf If \code{infPeriodSpec = "PS"} or \code{infPeriodSpec = "IDD"}, 
#' this will specify the maximum length of the infectious period.
#' @param dist If \code{infPeriodSpec = "PS"} this will specify the distribution 
#' for the path-specific removal probability. Currently only the exponential 
#' (\code{"exp"}), gamma (\code{"gamma"}), and Weibull (\code{"weibull"}) distributions 
#' are supported.
#' @param EKnown logical. If \code{EKnown = TRUE} the model will run assuming the
#' exposure times are known. \code{Estar} must be provided in \code{dat}.
#' @param WAIC logical. If \code{WAIC = TRUE}, WAIC will be computed using post burn-in samples.
#' @param progressBar logical. If \code{progressBar = TRUE} a progress bar will be displayed.
#' @param seed optional seeding for reproducibility.
#' 
#' @return A data frame of posterior samples or a two-element list with posterior samples and WAIC.
#' 
#' @details \code{mcmcSEIR} offers modelers three options to specify the infectious
#' period: the exponential distribution, the path-specific (PS) approach of Porter and Oleson (2013),
#' or the infectious-duration dependent (IDD) formulation of Ward et al. (upcoming).
#' 
#' All models will estimate a beta parameter associated with each column of the
#' matrix \code{X} and a rate parameter \code{rateE} associated with the exposure 
#' probability.
#' 
#' For the exponential model \code{infPeriodSpec = 'exp'}, only one additional parameter
#' will be estimated: \code{rateI} the rate parameter associated with the removal 
#' probability. Initial values and a prior for \code{rateI} need to be specified.
#' 
#' For the path-specific model \code{infPeriodSpec = 'PS'}, additional parameters 
#' corresponding to the distribution specified in \code{dist} will be estimated.
#' Initial values are named as \code{psParams}, which itself is a list with 
#' initial values for each associated parameter. Function giving priors for 
#' these parameters should be named \code{psParamsPrior} and 
#' contained in \code{priors}. For example, if \code{dist = 'gamma'}, then 
#' \code{inits$psParams = list(shape = , rate = )} will provide the initial 
#' values for the \code{shape} and \code{rate} parameters associated 
#' with the \code{gamma} distribution. \code{priors$psParamsPriors} should 
#' contain a function which computes the joint prior for these parameters.
#' 
#' For the infectious-duration dependent model \code{infPeriodSpec = 'IDD'}, 
#' additional parameters corresponding to the IDD curve specified in \code{iddFun}
#' will also be estimated. Initial values are named as \code{iddParams}, which 
#' itself is a list with initial values for each associated parameter. Function 
#' giving priors for these parameters should be named \code{iddParamsPrior} and 
#' contained in \code{priors}. For example, if \code{iddFun = logitIDD}, then
#' \code{inits$iddParams = list(rate = ..., mid = ...)} will provide the initial 
#' values and \code{priors$iddParamsPriors} should contain a function which computes
#' the joint prior for these parameters (see example).
#' 
#' If \code{WAIC = TRUE}, the WAIC will be calculated. WAIC is computed based on the 
#' log-likelihood values at each post burn-in iteration using the pWAIC2 
#' formulation described in Gelman et al. (2013), which multiplies the WAIC of 
#' Watanabe (2010) by 2 to be on the deviance scale.
#' 
#' 
#' @examples
#' \dontrun{
#' 
#' dat <- simSEIR(S0=999, E0=0, I0=1, N=1000, tau=100,
#'                beta=0.7, X=cbind(rep(1, 100)), rateE=0.1,
#'                infPeriodSpec = 'IDD', 
#'                infIDDParams = list(maxInf = 15,
#'                                    iddFun = logitIDD,
#'                                    iddParams = list(mid = 8, rate = 1.5)))
#' 
#' 
#' initsList <- list(beta = 1, rateE = 0.2, 
#'                   iddParams = list(mid = 4, rate = 0.1))
#'                   
#' priorsList <- list(betaPrior = function(x) dnorm(x, 0, 4, log = T),
#'                    rateEPrior = function(x) dgamma(x, 10, 100, log = T),
#'                    iddParamsPrior = function(x) {
#'                        dnorm(x["mid"], 5, 1, log = T) +
#'                        dgamma(x["rate"], 1, 1, log = T)
#'                    })
#' 
#' postRes <- mcmcSEIR(dat = dat, X=cbind(rep(1, 100)), 
#'                     inits = initsList, 
#'                     niter=10000, nburn=5000,
#'                     infPeriodSpec = 'IDD',
#'                     priors=priorsList,
#'                     iddFun = logitIDD, maxInf = 10,
#'                     progressBar = T)
#' 
#' }
#' @references 
#' Porter, Aaron T., and Oleson, Jacob J. "A path-specific SEIR model
#'  for use with general latent and infectious time distributions." \emph{Biometrics}
#'   69.1 (2013): 101-108.
#'   
#'   Gelman, Andrew, Hwang, Jessica and Vehtari, Aki. "Understanding predictive 
#'   information criteria for Bayesian models." \emph{Statistics and computing} 
#'   24.6 (2014): 997-1016.
#'   
#'   Watanabe, Sumio, and Opper, Manfred. "Asymptotic equivalence of Bayes cross 
#'   validation and widely applicable information criterion in singular learning
#'    theory." \emph{Journal of machine learning research} 11.12 (2010).
#' @export
mcmcSEIR <- function(dat, X, inits, niter, nburn,
                     infPeriodSpec = c('exp', 'PS', 'IDD'),
                     priors, iddFun = NULL, maxInf = NULL, dist = NULL,
                     EKnown = FALSE, WAIC = FALSE, progressBar = FALSE, 
                     seed = NULL) {
  
  infPeriodSpec <- match.arg(infPeriodSpec, c('exp', 'PS', 'IDD'))
  
  # check for valid MCMC specifications
  if (niter < 2) stop("Must run for at least 2 iterations")
  if (nburn > niter | nburn < 0) stop("invalid burn-in period")
  
  
  # check that data is valid
  if(!is.list(dat)) {
    stop('dat must be a list')
  }
  
  Istar <- dat$Istar
  S0 <- dat$S0
  E0 <- dat$E0
  I0 <- dat$I0
  N <- dat$N
  
  if (!EKnown) {
    if (!all(c('Istar', 'S0', 'E0', 'I0', 'N') %in% names(dat))) {
      stop('dat must contain Istar, S0, E0, I0, and N')
    }
  } else {
    Estar <- dat$Estar
    
    if (!all(c('Estar', 'Istar', 'S0', 'E0', 'I0', 'N') %in% names(dat))) {
      stop('dat must contain Estar, Istar, S0, E0, I0, and N')
    }
    
  }
  
  # check for valid seed and set if not specified
  if (is.null(seed)) {
    seed <- 1
  } else if (!is.numeric(seed)) {
    stop('seed must be numeric')
  }
  
  
  
  # check that X is valid
  # must be matrix (later check that beta inits are same length as ncol(X))
  if(!is.matrix(X)) {
    stop('X must be a matrix')
  }
  
  # must have one row for each time point
  if (nrow(X) != length(Istar)) {
    stop('X must have rows for each time point')
  }
  
  if (ncol(X) != length(inits$beta)) {
    stop('initial values for beta must be correct length')
  }
  
  # check that WAIC is a logical
  if (!is.logical(WAIC)) stop('WAIC must be logical value')
  
  if (infPeriodSpec == 'exp') {
    
    # check priors
    if(!is.list(priors)) {
      stop('priors must be a list')
    }
    
    if (!all(names(priors) %in% c('betaPrior', 'rateEPrior', 'rateIPrior'))) {
      stop('priors must be specified for beta, rateE, and rateI')
    }
    
    # extract priors from list
    betaPrior <- priors$betaPrior
    rateEPrior <- priors$rateEPrior
    rateIPrior <- priors$rateIPrior
    
    # check initial values
    if(!is.list(inits)) {
      stop('inits must be a list')
    }
    
    if (!all(names(inits) %in% c('beta', 'rateE', 'rateI'))) {
      stop('inits must be specified for beta, rateE, and rateI')
    }
    
    if (inits$rateE < 0 | inits$rateI < 0) {
      stop('invalid initial values')
    }
    
    betaProVar <- inits$beta/3
    
    maxInf <- ceiling(qexp(0.99, inits$rateI)) + 10
    
    
    if (!EKnown) {
      mcmcExp(Istar, X, maxInf, 
              S0, E0, I0, N, 
              niter, nburn, inits,
              betaProVar,
              betaPrior, rateEPrior, rateIPrior,
              WAIC, progressBar, seed)
      
    } else {
      mcmcExp_EKnown(Estar, Istar, X, maxInf, 
                     S0, E0, I0, N, 
                     niter, nburn, inits,
                     betaProVar,
                     betaPrior, rateEPrior, rateIPrior,
                     WAIC, progressBar, seed)
      
    }
    
    
  } else if (infPeriodSpec == 'PS') {
    
    # must specify dist and maxInf
    if (is.null(dist)) stop("For PS infectious period, dist must be specified")
    dist <- match.arg(dist, c('exp', 'gamma', 'weibull'))
    
    if (is.null(maxInf)) stop("For PS infectious period, maxInf must be specified")
    if (maxInf < 1) stop('maxInf must be >= 1')
    
    # check priors
    if(!is.list(priors)) {
      stop('priors must be a list')
    }
    
    if (!all(names(priors) %in% c('betaPrior', 'rateEPrior', 'psParamsPrior'))) {
      stop('priors must be specified for beta, rateE, and psParams')
    }
    
    # extract priors from list
    betaPrior <- priors$betaPrior
    rateEPrior <- priors$rateEPrior
    psParamsPrior <- priors$psParamsPrior
    
    # check initial values
    if(!is.list(inits)) {
      stop('inits must be a list')
    }
    
    if (!all(c('beta', 'rateE', 'psParams') %in% names(inits))) {
      stop('inits must be specified for beta, rateE, and psParams')
    }
    
    if (inits$rateE < 0) {
      stop('invalid initial values')
    }
    
    # initial values for psParams must match with specified distribution
    inits$psParams <- unlist(psParamsCheck(inits$psParams, dist) )
    
    betaProVar <- inits$beta/3
    psParamsProVar <- unlist(inits$psParams) / 4
    
    if (!EKnown) {
      mcmcPS(Istar, X, maxInf, S0, E0, I0, N, 
             niter, nburn, inits,
             dist, 
             betaProVar, psParamsProVar,
             betaPrior, rateEPrior, psParamsPrior,
             WAIC, progressBar, seed)
    } else {
      mcmcPS_EKnown(Estar, Istar, X, maxInf, S0, E0, I0, N, 
                    niter, nburn, inits,
                    dist, 
                    betaProVar, psParamsProVar,
                    betaPrior, rateEPrior, psParamsPrior,
                    WAIC, progressBar, seed)
    }
    
    
  } else if (infPeriodSpec == 'IDD') {
    
    # iddFun must be specified
    if (is.null(iddFun)) stop("For IDD infectious period, iddFun must be specified")
    
    # check initial values
    if(!is.list(inits)) {
      stop('inits must be a list')
    }
    
    if (!all(c('beta', 'rateE', 'iddParams') %in% names(inits))) {
      stop('inits must be specified for beta, rateE, and iddParams')
    }
    
    if (inits$rateE < 0) {
      stop('invalid initial values')
    }
    
    
    # if spline model, need to fix the XBasis argument and remove it from the parameters
    # iddFun could be character or function
    if (is.character(iddFun) & is.logical(all.equal(get(iddFun), splineIDD))) {
      XBasis <- inits$iddParams$XBasis
      inits$iddParams <- inits$iddParams[-which(names(inits$iddParams) == 'XBasis')]
      iddFun <- fixFunArgs(substitute(splineIDD(XBasis=XBasis)))
    } else if (is.function(iddFun) & is.logical(all.equal(iddFun, splineIDD))) {
      XBasis <- inits$iddParams$XBasis
      inits$iddParams <- inits$iddParams[-which(names(inits$iddParams) == 'XBasis')]
      iddFun <- fixFunArgs(substitute(splineIDD(XBasis=XBasis)))
    }
    
    iddParamNames <- do.call(iddFun, args = list(x = 1,
                                                 params = inits$iddParams,
                                                 getParams = T))
    
    tryCatch(do.call(iddFun, args = list(x = 1:maxInf,
                                         params = inits$iddParams)),
             error = function(e) stop('invalid initial values'))
    
    # check maxInf 
    if (is.null(maxInf)) stop("For IDD infectious period, maxInf must be specified")
    if (maxInf < 1) stop('maxInf must be >= 1')
    
    # check priors
    if(!is.list(priors)) {
      stop('priors must be a list')
    }
    
    if (!all(names(priors) %in% c('betaPrior', 'rateEPrior', 'iddParamsPrior'))) {
      stop('priors must be specified for beta, rateE, and iddParams')
    }
    
    # extract priors from list
    betaPrior <- priors$betaPrior
    rateEPrior <- priors$rateEPrior
    iddParamsPrior <- priors$iddParamsPrior
    
    
    
    betaProVar <- inits$beta/3
    iddParamsProVar <- unlist(inits$iddParams)/3
    
    if (!EKnown) {
      mcmcIDD(Istar, X, maxInf, S0, E0, I0, N, 
              niter, nburn, inits,
              iddFun, iddParamNames,
              betaProVar, iddParamsProVar,
              betaPrior, rateEPrior, iddParamsPrior,
              WAIC, progressBar, seed)
    } else {
      mcmcIDD_EKnown(Estar, Istar, X, maxInf, S0, E0, I0, N, 
                     niter, nburn, inits,
                     iddFun, iddParamNames,
                     betaProVar, iddParamsProVar,
                     betaPrior, rateEPrior, iddParamsPrior,
                     WAIC, progressBar, seed)
    }
    
    
    
    
  }
  
  
  
  
}
