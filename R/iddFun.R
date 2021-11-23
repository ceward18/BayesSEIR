################################################################################
# IDD functions
# must only have 3 parameters - x, params, and getParams (used internally)
################################################################################

#' IDD transmissibility curve using the gamma PDF
#'
#' @param x A vector ranging from 1 to the maximum length of the infectious period.
#' @param params A list giving parameter values for the shape and rate of the gamma PDF.
#' @param getParams logical. Tells the function to return the parameter names.
#' @return Gamma PDF  
#' \code{dgamma(x, shape = params$shape, rate = params$rate)}
#' @examples
#' dgammaIDD(1:15, params = list(shape = 4, rate = 1))
#' @export
dgammaIDD <- function(x, params, getParams = F) {
  
  if (getParams) {
    return(c('shape', 'rate'))
  }
  
  if (!all(names(params) %in% c('shape', 'rate'))) {
    stop('params must contain shape and rate')
  }
  
  shape <- params$shape
  rate <- params$rate
  
  # check that shape and rate are positive 
  # this will cause negative proposals to be automatically rejected
  if (shape < 0) stop("shape must be positive")
  if (rate < 0) stop("rate must be positive")
  
  dgamma(x, shape, rate)
}


#' IDD transmissibility curve using the log normal PDF
#'
#' @param x A vector ranging from 1 to the maximum length of the infectious period.
#' @param params A list giving parameter values for the meanlog and sdlog of the log normal PDF.
#' @param getParams logical. Tells the function to return the parameter names.
#' @return Log normal PDF 
#' \code{dlnorm(x, meanlog = params$meanlog, sdlog = params$sdlog)}
#' @examples
#' dlnormIDD(1:15, params = list(meanlog = 2, sdlog = 1))
#' @export
dlnormIDD <- function(x, params, getParams = F) {
  
  if (getParams) {
    return(c('meanlog', 'sdlog'))
  }
  
  if (!all(names(params) %in% c('meanlog', 'sdlog'))) {
    stop('params must contain meanlog and sdlog')
  }
  
  meanlog <- params$meanlog
  sdlog <- params$sdlog
  
  if (sdlog < 0) stop("sdlog must be positive")
  
  dlnorm(x, meanlog, sdlog)
}

#' IDD transmissibility curve using logistic decay
#'
#' @param x A vector ranging from 1 to the maximum length of the infectious period.
#' @param params A list giving parameter values for the decay rate \code{rate}
#'  and inflection point \code{mid} of the logistic decay function.
#' @param getParams logical. Tells the function to return the parameter names.
#' @return \code{1 / (1 + exp(params$rate * (x - params$mid)))}
#' @examples
#' logitIDD(1:15, params = list(mid = 8, rate = 1.5))
#' @export
logitIDD <- function(x, params, getParams = F) {
  
  if (getParams) {
    return(c('rate', 'mid'))
  }
  
  if (!all(names(params) %in% c('rate', 'mid'))) {
    stop('params must contain rate and mid')
  }
  
  rate <- params$rate
  mid <- params$mid
  
  1 / (1 + exp(rate * (x - mid)))
}

#' IDD transmissibility curve using a spline basis
#'
#' @param x A vector ranging from 1 to the maximum length of the infectious period.
#' @param params A list giving coefficients for each basis in the basis matrix, 
#' named as \code{b1}, \code{b2}, ...
#' @param XBasis Basis matrix. Generally the output from something such as \code{splines::bs()}.
#' @param getParams logical. Tells the function to return the parameter names.
#' @return \code{exp(XBasis %*% b)}, where \code{b = c(params$b1, params$b2, ...)} 
#' 
#' @details The spline IDD curve is exponentiated to ensure f(x) > 0
#' 
#' 
#' @examples
#' splineIDD(1:15, params = list(b1 = 3, b2 = 1, b3 = -3),
#'                 XBasis = splines::bs(1:15, df = 3))
#' @export
splineIDD <- function(x, params, XBasis, getParams = F) {
  
  if (getParams) {
    return(paste0('b', 1:ncol(XBasis)))
  }
  
  if (!all(names(params) %in% c(paste0('b', 1:ncol(XBasis))))) {
    stop('params must contain XBasis and associated b coefficients')
  }
  
  # convert to vector
  b <- unlist(params)
  
  if (length(b) != ncol(XBasis)) {
    stop('params must contain b for each basis')
  }

  c(exp(XBasis %*% b))
}








