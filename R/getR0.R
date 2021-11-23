################################################################################
# Functions to compute basic reproductive numbers
################################################################################


#' Function to compute the basic R0 given a set of parameters.
#'
#' @param infPeriodSpec A character indicating which model for the infectious period
#' will be used. 
#' @param beta vector with elements for each column of \code{X} which specify the 
#' transmission rate over time.
#' @param X The design matrix corresponding to the intensity process over epidemic time.
#' @param N Integer, population size.
#' @param infExpParams A \emph{list} giving parameters specific to the exponentially
#' distributed infectious period. See details for specifics.
#' @param infPSParams A \emph{list} giving parameters specific to the path-specific
#' distributed infectious period. See details for specifics.
#' @param infIDDParams A \emph{list} giving parameters specific to the infectious
#' duration-dependent infectious period. See details for specifics.
#' 
#' @return A vector with the value of the basic reproductive number R0(t) over the
#' epidemic
#' 
#' @details \code{getR0} computes the basic reproductive number which may change
#' over epidemic time due to an intervention. Reproductive numbers can be calculated
#' under three specifications of the infectious period: the exponential distribution, 
#' the path-specific (PS) approach of Porter and Oleson (2013),
#' or the infectious-duration dependent (IDD) formulation of Ward et al. (upcoming).
#' 
#' \code{infPeriodSpec} determines which model will be used to simulate the epidemic.
#' Model specific parameter values are entered using either \code{infExpParams}, 
#' \code{infPSParams}, or \code{infIDDParams}.
#' 
#' All models use the \code{beta} parameter vector associated with each column of the
#' matrix \code{X} to describe transmission over epidemic time. All models 
#' assume exponentially distributed time spent in the latent period with rate 
#' parameter, \code{rateE}.
#' 
#' For the exponential model \code{infPeriodSpec = 'exp'}, only one additional 
#' parameter needs to be specified in \code{infExpParams}: \code{rateI}, which is
#' the rate parameter associated with the removal probability. 
#' 
#' For the path-specific model \code{infPeriodSpec = 'PS'}, \code{infPSParams}
#' should contain three elements: \code{dist}, \code{maxInf}, and \code{psParams}.
#' \code{dist} gives the distribution used to describe the length of time spent 
#' in the infectious period Currently only the exponential (\code{"exp"}), 
#' gamma (\code{"gamma"}), and Weibull (\code{"weibull"}) distributions are supported.
#' \code{maxInf} corresponds to the maximum length of time an individual can 
#' spend in the infectious compartment, at which point the removal probability 
#' becomes 1. \code{psParams} must provide a \emph{list} of parameter values for the parameters
#' associated with the chosen distribution. For \code{dist = 'gamma'}, these are 
#' the \code{shape} and \code{rate} parameters and for \code{dist = 'weibull'}, 
#' these are the \code{shape} and \code{scale} parameters.
#' 
#' For the infectious-duration dependent model \code{infPeriodSpec = 'IDD'}, 
#' \code{infIDDParams} should contain three elements: \code{maxInf}, \code{iddFun}, 
#' and \code{iddParams}. \code{maxInf} corresponds to the total length of time 
#' each individual spends in the infectious compartment. \code{iddFun} gives
#' the IDD function used to describe the IDD curve. \code{iddParams} must provide 
#' a \emph{list} of parameter values for the parameters associated with the chosen \code{iddFun}. 
#' For example, if \code{iddFun = dgammaIDD}, these are the \code{shape} and 
#' \code{rate} parameters.
#' 
#' @examples
#' 
#' # simulate epidemic over 100 days
#' tau <- 100
#' 
#' # specify the design matrix so there is a change point in transmission at time 50
#' X <- cbind(1, c(rep(0, 49), rep(1, tau - 49)))
#' 
#' 
#' # simulate using exponentially distributed infectious period
#' datExp <- simSEIR(S0 = 999, E0 = 0, I0 = 1, N = 1000, tau = tau,
#'                   beta = c(0.1, -2), X = X, rateE = 0.1,
#'                   infPeriodSpec = 'exp',
#'                   infExpParams = list(rateI = 0.2))
#'                   
#' # plot incidence curve
#' plot(datExp$Istar, type = 'l', main = 'Incidence', xlab = 'Epidemic Time',
#'      ylab = 'Count')
#'      
#' # get reproductive number over time
#' getR0(infPeriodSpec = 'exp', beta = c(0.1, -2), X = X, N = 1000,
#'       infExpParams = list(rateI = 0.2))
#'      
#' # simulate using IDD infectious period using the gamma PDF
#' datIDD <- simSEIR(S0 = 999, E0 = 0, I0 = 1, N = 1000, tau = tau,
#'                   beta = c(0.7, -1.5), X = X, rateE = 0.1,
#'                   infPeriodSpec = 'IDD',
#'                   infIDDParams = list(maxInf = 14,
#'                                       iddFun = dgammaIDD,
#'                                       iddParams = list(shape = 4, rate = 1))) 
#'                                       
#' # plot incidence curve
#' plot(datIDD$Istar, type = 'l', main = 'Incidence', xlab = 'Epidemic Time',
#'      ylab = 'Count')
#'      
#' # get reproductive number over time
#' getR0(infPeriodSpec = 'IDD', beta = c(0.7, -1.5), X = X, N = 1000,
#'       infIDDParams = list(maxInf = 14,
#'                           iddFun = dgammaIDD,
#'                           iddParams = list(shape = 4, rate = 1)))
#'
#' @references Porter, Aaron T., and Jacob J. Oleson. "A path-specific SEIR model
#'  for use with general latent and infectious time distributions." \emph{Biometrics}
#' @export
getR0 <- function(infPeriodSpec = c('exp', 'PS', 'IDD'), beta, X, N,
                  infExpParams = NULL,
                  infPSParams = NULL,
                  infIDDParams = NULL) {
  
  infPeriodSpec <- match.arg(infPeriodSpec, c('exp', 'PS', 'IDD'))
  
  # check N is valid
  if (N < 1) {
    stop("Population size, N, must be positive")
  }
  
  
  # check that beta and X are valid
  # must be matrix (later check that beta inits are same length as ncol(X))
  if(!is.matrix(X)) {
    stop('X must be a matrix')
  }
  
  if (ncol(X) != length(beta)) {
    stop('Value for beta must be correct length (one element for each column of X)')
  }
  
  
  if (infPeriodSpec == 'exp') {
    
    if (is.null(infExpParams)) stop('infExpParams must be specified')
    
    if (!is.list(infExpParams)) {
      stop('infExpParams must be a list')
    }
    
    if (!all(names(infExpParams) %in% c('rateI'))) {
      stop("if infPeriodSpec = 'exp', infExpParams must contain 'rateI'")
    }
    
    rateI <- infExpParams$rateI
    
    if (rateI < 0) stop("rateI must be positive")
    
    getR0Exp(X, beta, N, rateI) 
    
  } else if (infPeriodSpec == 'PS') {
    
    if (is.null(infPSParams)) stop('infPSParams must be specified')
    
    if (!is.list(infPSParams)) {
      stop('infPSParams must be a list')
    }
    
    if (!all(names(infPSParams) %in% c('dist', 'psParams', 'maxInf'))) {
      stop("if infPeriodSpec = 'PS', infPSParams must contain 'dist', 'psParams', 'maxInf'")
    }
    
    dist <- match.arg(infPSParams$dist, c('exp', 'gamma', 'weibull'))
    psParams <- infPSParams$psParams
    maxInf <- infPSParams$maxInf
    
    psParams <- psParamsCheck(psParams, dist) 

    if (maxInf < 1) stop('maxInf must be >= 1')
    
    getR0PS(X, beta, N, maxInf, dist, psParams) 
  
    
  } else if (infPeriodSpec== 'IDD') {
    
    if (is.null(infIDDParams)) stop('infIDDParams must be specified')
    
    if (!is.list(infIDDParams)) {
      stop('infIDDParams must be a list')
    }
    
    if (!all(names(infIDDParams) %in% c('iddFun', 'iddParams', 'maxInf'))) {
      stop("if infPeriodSpec = 'IDD', infIDDParams must contain 'iddFun', 'iddParams', 'maxInf'")
    }
    
    iddFun <- infIDDParams$iddFun
    iddParams <- infIDDParams$iddParams
    maxInf <- infIDDParams$maxInf
    

    # if spline model, need to fix the XBasis argument and remove it from the parameters
    if (!is.character(all.equal(iddFun, splineIDD))) {
      
      XBasis <- iddParams$XBasis
      iddParams <- iddParams[-which(names(iddParams) == 'XBasis')]
      iddFun <- fixFunArgs(substitute(splineIDD(XBasis=XBasis)))
    }
    
    if (maxInf < 1) stop('maxInf must be >= 1')
    
    IDDCurve <- do.call(iddFun, args = list(x = 1:maxInf,
                                            params = iddParams))
    
    
    getR0IDD(X, beta, N, IDDCurve, maxInf)
  } 
  
  
}


# function to get R0 for exponential model
getR0Exp <- function(X, beta, N, rateI) {
  
  maxInf <- 100
  
  # extend X so R0 can be calculated for the last maxInf time points
  X_ext <- rbind(X, matrix(rep(X[nrow(X),], maxInf + 1), ncol = ncol(X), byrow = T))
  
  # probability of transition given 1 infectious individual
  pi_SE <- c(1 - exp(- exp(X_ext %*% beta) / N ))
  
  # infectious probability of removal
  pi_IR <- 1 - exp(-rateI)
  
  # for exponential periods
  stayProb <- rep(1 - pi_IR, maxInf)
  multVec <- cumprod(stayProb)
  multVec <- c(1, multVec)
  
  addUpSmooth(pi_SE, multVec, N)
}


# function to get R0 for PS model
getR0PS <- function(X, beta, N, maxInf, dist, psParams) {
  
  # extend X so R0 can be calculated for the last maxInf time points
  X_ext <- rbind(X, matrix(rep(X[nrow(X),], maxInf + 1), ncol = ncol(X), byrow = T))
  
  # probability of transition given 1 infectious individual
  pi_SE <- c(1 - exp(- exp(X_ext %*% beta) / N ))
  
  # infectious probability of removal
  pi_IR <- psProbVec(maxInf, dist, psParams) 
  stayProb <- 1 - pi_IR
  multVec <- cumprod(stayProb)
  multVec <- c(1, multVec)
  
  addUpSmooth(pi_SE, multVec, N)
}

# function to get R0 for IDD model
getR0IDD <- function(X, beta, N, IDDCurve, maxInf) {
  
  # extend X so R0 can be calculated for the last maxInf time points
  X_ext <- rbind(X, matrix(rep(X[nrow(X),], maxInf), ncol = ncol(X), byrow = T))
  
  theta <-  c(exp(X_ext %*% beta) )
  
  addUpSmoothIDD(theta, IDDCurve, N, maxInf)
}


# helper function to sum forward in time for IDD
addUpSmoothIDD <- function(theta, IDDCurve, N, maxInf){
  nTime <- length(theta)
  bw <- maxInf
  sumSmooth <- rep(NA, nTime - bw)
  for(k in 1:(nTime - bw)){
    t1 <- k
    t2 <- k + bw - 1
    pi_SE <- 1 - exp(-theta[t1:t2] * IDDCurve / N)
    
    sumSmooth[k] <- sum(pi_SE) * N
  }
  sumSmooth
}

# helper function to sum forward in time for exp or PS
addUpSmooth <- function(pi_SE, multVec, N){
  nTime <- length(pi_SE)
  bw <- length(multVec)
  sumSmooth <- rep(NA, nTime - bw)
  for(k in 1:(nTime - bw)){
    t1 <- k
    t2 <- k + bw - 1
    sumSmooth[k] <- sum(pi_SE[t1:t2] * multVec) * N
  }
  sumSmooth
}