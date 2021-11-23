
################################################################################
# compiled functions for likelihoods
################################################################################


# calculate probability of transmission at each time/location
# exponential and PS probability
p_SE <- function(X, beta, infVec, popSizesInv) {
  
  theta <- exp(X %*% beta) * infVec * popSizesInv 
  
  # compute probability
  1 - exp(-theta - 1e-10)
}


# IDD transmission probability
p_SE_IDD <- function(X, beta, popSizesInv, iddFun, iddParams, fullI, maxInf) {
  
  # theta for each time and location
  theta_components <- exp(X %*% beta)
  
  # IDD transmissibility curve
  iddCurve <- do.call(iddFun, args = list(x = 1:maxInf,
                                            params = as.list(iddParams)))
  
  do.call(dgamma, args = list(x = 1:20, shape = 4, rate = 1))
  
  delta <- colByVec(fullI, iddCurve) * popSizesInv 
  
  theta <- theta_components * delta 
  
  # compute probability
  1 - exp(-theta - 1e-10)
}


# binomial draw for transmission with constant transmissibility
f_SE <- function(susVec, Estar, X, beta, infVec, popSizesInv) {
  
  pi_SE <- p_SE(X=X, beta=beta, infVec=infVec, popSizesInv=popSizesInv)
  
  sum(dbinom(Estar, susVec, c(pi_SE), log = T))
  
}


# binomial draw for transmission with IDD transmissibility
f_SE_IDD <- function(susVec, Estar, X, beta, popSizesInv, 
                     iddFun, iddParams, fullI, maxInf) {
  
  
  pi_SE <- p_SE_IDD(X=X, beta=beta,popSizesInv=popSizesInv, 
                        iddFun=iddFun, iddParams=iddParams,
                        fullI=fullI, maxInf=maxInf)
  
  sum(dbinom(Estar, susVec, c(pi_SE), log = T))
  
}

# binomial draw for exponential transitions with rate parameter
# E to I compVec = expVec, transVec = Istar
# I to R compVec = infVec, transVec = Rstar
f_exp <- function(compVec, transVec, rate) {
  probT <- 1 - exp(-rate)
  sum(dbinom(transVec, compVec, probT, log = T))
}

# path-specific transition from I to R
# can use exponential, gamma or weibull distributions
f_PS <- function(fullMat, fullStar, dist, psParams, maxLength) {
  # compute possible transition probabilities based on time infected
  
  pi_IRVec <- psProbVec(maxInf=maxLength, dist=dist, psParams=as.list(psParams))
  pi_IRVec[zapsmall(pi_IRVec)==1] <- 0.9999
  
  # sum over locations and time
  fullVec <- colSums(fullMat)
  starVec <- colSums(fullStar)
  
  maxTimeObs <- max(length(fullVec), length(starVec))
  
  sum(lchoose(fullMat[,1:ncol(fullMat)], fullStar[,1:ncol(fullStar)])) + 
    sum(starVec*log(pi_IRVec[1:maxTimeObs]) + (fullVec-starVec)*log(1 - pi_IRVec[1:maxTimeObs]))

}


################################################################################
# full log-likelihoods for WAIC calculations

# exponential model
ll_exp <- function(susVec, Estar, X, beta, infVec, popSizesInv, 
               expVec, Istar, rateE, Rstar, rateI) {
  
  pi_SE <- p_SE(X=X, beta=beta, infVec=infVec, popSizesInv=popSizesInv)
  
  pi_EI <- 1 - exp(-rateE)
  pi_IR <- 1 - exp(-rateI)
  
  SE_lik <- dbinom(Estar, susVec, c(pi_SE), log = T)
  EI_lik <- dbinom(Istar, expVec, pi_EI, log = T)
  IR_lik <- dbinom(Rstar, infVec, pi_IR, log = T)
  
  SE_lik + EI_lik + IR_lik
  
}

# path-specific model
ll_PS <- function(susVec, Estar, X, beta, infVec, popSizesInv, 
                   expVec, Istar, rateE, maxInf, dist, psParams, 
                  fullI, fullRstar) {
  
  pi_SE <- p_SE(X=X, beta=beta, infVec=infVec, popSizesInv=popSizesInv)
  
  pi_EI <- 1 - exp(-rateE)
  
  pi_IRVec <- psProbVec(maxInf=maxInf, dist=dist, psParams=as.list(psParams)) 
  pi_IRVec[zapsmall(pi_IRVec)==1] <- 0.9999
  
  SE_lik <- dbinom(Estar, susVec, c(pi_SE), log = T)
  EI_lik <- dbinom(Istar, expVec, pi_EI, log = T)
  
  IR_lik <- rep(0, length(EI_lik))
  for (j in 1:ncol(fullI)) {
    IR_lik <- IR_lik + dbinom(fullRstar[,j], fullI[,j], pi_IRVec[j], log = T)
  }
  
  SE_lik + EI_lik + IR_lik
  
}

# IDD model
ll_IDD <- function(susVec, Estar, X, beta, popSizesInv, iddFun, iddParams,
                   fullI, maxInf, expVec, Istar, rateE) {
  
  pi_SE <- p_SE_IDD(X=X, beta=beta, popSizesInv=popSizesInv, 
                    iddFun=iddFun, iddParams=iddParams,
                    fullI=fullI, maxInf=maxInf)
  
  pi_EI <- 1 - exp(-rateE)
  
  SE_lik <- dbinom(Estar, susVec, c(pi_SE), log = T)
  EI_lik <- dbinom(Istar, expVec, pi_EI, log = T)
  
  SE_lik + EI_lik
  
}


