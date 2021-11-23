################################################################################
# simulations scripts (called by simSEIR)
# simExp - exponential infectious period
# simPS - path-specific infectious period
# simIDD - IDD infectious period
################################################################################

# exponential infectious period
simExp <- function(S0, E0, I0, N, tau,
                   beta, X, rateE, rateI) {
  
  # initialization of compartments
  S <- E <- I <- R <- Estar <- Istar <- Rstar <- rep(NA, tau + 1)
  S[1] <- S0 
  E[1] <- E0
  I[1] <- I0
  R[1] <- N - S[1] - E[1] - I[1]
  
  # initialize pi_EI and pi_IR
  pi_EI <- 1 - exp(-rateE)
  pi_IR <- 1 - exp(-rateI)
  
  # transmission may change over epidemic time
  theta_components <- exp(X %*% beta)
  
  for (i in 1:tau){
    pi_SE <- 1 - exp(- theta_components[i] * I[i] / N)
    
    Estar[i] <- rbinom(n = 1, size = S[i], prob = pi_SE)
    Istar[i] <- rbinom(n = 1, size = E[i], prob = pi_EI)
    Rstar[i] <- rbinom(n = 1, size = I[i], prob = pi_IR)
    
    # update compartments according to transitions
    S[i + 1] <- S[i] - Estar[i]
    E[i + 1] <- E[i] + Estar[i] - Istar[i]
    I[i + 1] <- I[i] + Istar[i] - Rstar[i]
    R[i + 1] <- R[i] + Rstar[i]
    
  }
  
  list(N = N, S0 = S0, E0 = E0, I0 = I0,
       S = S[-(tau + 1)], E = E[-(tau + 1)], I = I[-(tau + 1)], R = R[-(tau + 1)],
       Estar = Estar[-(tau + 1)], Istar = Istar[-(tau + 1)], Rstar = Rstar[-(tau + 1)])
  
}

# path-specific infectious period
# possible distributions are exponential, gamma, or Weibull
simPS <- function(S0, E0, I0, N, tau,
                   beta, X, rateE, 
                  dist, psParams, maxInf) {
  
  # initialization of compartments
  S <- R <- E <- Estar <- Istar <- rep(NA, tau + 1)
  S[1] = S0 
  E[1] = E0
  I <- Rstar <- matrix(0, nrow = tau + 1, ncol = maxInf)
  I[1,1] = I0
  R[1] = N - S[1] - E[1] - I[1,1]
  
  # initialize pi_EI and pi_IR
  pi_EI <- 1 - exp(-rateE)
  pi_IRVec <- psProbVec(maxInf, dist, psParams) 

  # transmission may change over epidemic time
  theta_components <- exp(X %*% beta)
  
  for (i in 1:tau){
    
    pi_SE <- 1 - exp(- theta_components[i] * sum(I[i,]) / N)
    
    Estar[i] <- rbinom(n = 1, size = S[i], prob = pi_SE)
    Istar[i] <- rbinom(n = 1, size = E[i], prob = pi_EI)
    Rstar[i,] <- rbinom(n = rep(1, maxInf), size = I[i,], prob = pi_IRVec)
    
    # update compartments according to transitions
    S[i + 1] <- S[i] - Estar[i]
    E[i + 1] <- E[i] + Estar[i] - Istar[i]
    I[i + 1,] <- c(0, I[i,] - Rstar[i,])[-(maxInf + 1)]
    I[i + 1,1] <- I[i + 1,1] + Istar[i]
    R[i + 1] <- R[i] + sum(Rstar[i,])
    
  }
  
  list(N = N, S0 = S0, E0 = E0, I0 = I0,
       S = S[-(tau + 1)], E = E[-(tau + 1)], I = I[-(tau + 1),], R = R[-(tau + 1)],
       Estar = Estar[-(tau + 1)], Istar = Istar[-(tau + 1)], Rstar = Rstar[-(tau + 1),])
  
}

# IDD infectious period
simIDD <- function(S0, E0, I0, N, tau,
                  beta, X, rateE, 
                  IDDCurve, maxInf) {
  
  # initialization of compartments
  S <- E <- R <- Estar <- Istar <- Rstar <- rep(0, tau + 1)
  S[1] = S0 
  E[1] = E0
  I <- matrix(0, nrow = tau + 1, ncol = maxInf)
  I[1,1] = I0
  R[1] = N - S[1] - E[1] - I[1,1]
  Rstar[1 + maxInf] <- I0
  
  # initialize pi_EI and pi_IR
  pi_EI <- 1 - exp(-rateE)
  
  # transmission may change over epidemic time
  theta_components <- exp(X %*% beta)
  
  for (i in 1:tau){
    
    pi_SE <- 1 - exp(- theta_components[i] * sum(I[i,] * IDDCurve) / N)
    
    Estar[i] <- rbinom(n = 1, size = S[i], prob = pi_SE)
    Istar[i] <- rbinom(n = 1, size = E[i], prob = pi_EI)
    
    # Rstar is fully determined by Istar (fixed length infectious period)
    if ((i + maxInf - 1) <= tau) {
      Rstar[i + maxInf] <- Istar[i]
    }
    
    # update compartments according to transitions
    S[i + 1] <- S[i] - Estar[i]
    E[i + 1] <- E[i] + Estar[i] - Istar[i]
    I[i + 1,] <- c(0, I[i,])[-(maxInf + 1)]
    I[i + 1, 1] <- Istar[i]
    R[i + 1] <- R[i] + Rstar[i]
    
  }
  
  list(N = N, S0 = S0, E0 = E0, I0 = I0,
       S = S[-(tau + 1)], E = E[-(tau + 1)], I = I[-(tau + 1),], R = R[-(tau + 1)],
       Estar = Estar[-(tau + 1)], Istar = Istar[-(tau + 1)], Rstar = Rstar[-(tau + 1)])
  
}
