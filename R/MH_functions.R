################################################################################
# Metropolis Hastings functions
# adaptive and non-adaptive MH functions
# proposal distribution functions
################################################################################

# single metropolis update (adaptive)
# Vihola 2012
MH_update_adapt <- function(x0, f, S, currentIter, nBurn, ...) {
  
  p <- length(x0)
  
  U <- rnorm(p)
  x1 <- round(c(x0 + S %*% U), 6)
  
  names(x1) <- names(x0)
  
  # error handling for proposals outside the valid parameter domain
  fcOld <- f(x = x0, ...)
  
  # this will cause exp(a) = 0, and u is always > -Inf, i.e. proposal will always
  #   be rejected
  fcNew <- tryCatch(f(x = x1, ...), error = function(e) -Inf)
  
  a <- sum(fcNew - fcOld)
  u <- log(runif(1, 0, 1))
  
  # accept or reject proposal
  if (u < a) x <- x1
  else x <- x0
  
  if (currentIter < nBurn) {
    S <- ramcmc::adapt_S(S, U, min(1, exp(a)), currentIter - 1)
  }
  
  list(x = x, S = S)
  
}


# single metropolis hastings update 
MH_update <- function(x0, f, g, g_s, sd, ...){
  # Multivariate version of the MH
  x1 <- g_s(x0, sd=sd)
  
  a <- sum((f(x=x1,...) - f(x=x0,...))) + 
    g(x1, x0, sd=sd) - g(x0, x1, sd=sd)
  u <- log(runif(1,0,1))
  
  
  if (u < a){ 
    return(x1)
  }
  return(x0)
}

# truncated normal proposal (can't be negative)
g_s_positive <- function(from, sd) {
  truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = from, sd = sd)
}
g_positive <- function(from, to, sd) {
  sum(log(truncnorm::dtruncnorm(to, a = 0, b = Inf, mean=from, sd = sd)))
}

# continuous proposal functions
g_s <- function(from, sd) {
  MASS::mvrnorm(1, mu=from, Sigma = diag(sd, length(from)))
}
g <- function(from, to, sd) {
  mvtnorm::dmvnorm(to, mean=from, sigma = diag(sd, length(from)), log = T)
}



