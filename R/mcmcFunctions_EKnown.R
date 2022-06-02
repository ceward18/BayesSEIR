################################################################################
# MCMC scripts for known exposure period
# initialize and run chains for specified number of iterations
# mcmcExp_EKnown
# mcmcPS_EKnown
# mcmcIDD_EKnown
################################################################################


################################################################################
### mcmcExp
################################################################################
mcmcExp_EKnown <- function(Estar, Istar, X, maxInf, S0, E0, I0, N, 
                    niter, nburn, inits,
                    betaProVar,
                    betaPrior, rateEPrior, rateIPrior,
                    WAIC, progress) {
  
  # useful things from data
  epiSize <- sum(Istar)
  popSizesInv <- N^(-1)
  
  Rstar <- c(I0, Istar[-length(Istar)])
  
  # initalize compartment vectors
  susVec <- getS(S0, Estar)
  expVec <- getE(E0, Estar, Istar)
  infVec <- getI(I0, Istar, Rstar)
  
  # initialize storage
  # variables in transmission probability
  transProbPost <- matrix(NA, nrow = niter, ncol = ncol(X))
  rateEPost <- rep(NA, niter)
  rateIPost <- rep(NA, niter)
  
  transProbPost[1,] <- inits$beta
  rateEPost[1] <- inits$rateE
  rateIPost[1] <- inits$rateI
  
  # initialize covariance for transmission probability proposal
  S <- diag(betaProVar, ncol(transProbPost))
  
  # initialize storage of log-likelihood if WAIC is requested
  if (WAIC) {
    ll_store <- matrix(NA, nrow = length(Istar), ncol = niter)
    ll_store[,1] <- ll_exp(susVec=susVec, Estar=Estar, X=X, beta=transProbPost[1,],
                           infVec=infVec, popSizesInv=popSizesInv, 
                           expVec=expVec, Istar=Istar, rateE=rateEPost[1],
                           Rstar=Rstar, rateI=rateIPost[1]) 
  }
  
  # create progress bar
  if (progress) {
    pb <- txtProgressBar(min = 1, max = niter, style = 3)
  }
  
  for (i in 2:niter) {
    
    # transmission probability
    transProbUpdate <- MH_update_adapt(transProbPost[i - 1, ], fcTransProb, S, i, nburn,
                                       susVec=susVec, Estar=Estar, X=X,  
                                       infVec=infVec, popSizesInv=popSizesInv,
                                       betaPrior=betaPrior)
    
    transProbPost[i,] <- transProbUpdate$x
    S <- transProbUpdate$S
    
    # exponential rate parameter for E to I
    rateEPost[i] <- MH_update(rateEPost[i-1], fcRate, g_positive, g_s_positive, sd=0.05, 
                              compVec=expVec, transVec=Istar, ratePrior=rateEPrior)
    
    # exponential rate parameter for I to R
    rateIPost[i] <- MH_update(rateIPost[i-1], fcRate, g_positive, g_s_positive, sd=0.05, 
                              compVec=infVec, transVec=Rstar, ratePrior=rateIPrior)
    
   
    # Rstar
    RstarProp <- update_Rstar_exp(Rstar=Rstar, maxInf=maxInf, epiSize=epiSize,
                                  susVec=susVec, Estar=Estar, X=X,
                                  beta=transProbPost[i,], infVec=infVec, 
                                  popSizesInv=popSizesInv, rateI=rateIPost[i], 
                                  expVec=expVec, I0=I0, Istar=Istar, percent=0.1) 
    
    
    Rstar <- RstarProp$Rstar
    infVec <- RstarProp$infVec
    
    # update log-likelihood
    if (WAIC) {
      ll_store[,i] <- ll_exp(susVec=susVec, Estar=Estar, X=X, beta=transProbPost[i,],
                             infVec=infVec, popSizesInv=popSizesInv, 
                             expVec=expVec, Istar=Istar, rateE=rateEPost[i],
                             Rstar=Rstar, rateI=rateIPost[i]) 
      
    }
    
    # update progress bar
    if (progress) {
      setTxtProgressBar(pb, i)
    }
    
  }
  
  # export posteriors in one data frame
  fullPost <- cbind.data.frame(transProbPost, rateEPost, rateIPost)
  colnames(fullPost) <- c(paste0('beta', 1:ncol(transProbPost)), 'rateE', 'rateI')
  
  if (WAIC) {
    return(list(fullPost = fullPost[(nburn + 1):niter,], 
                WAIC = getWAIC(ll_store[,(nburn + 1):niter])))
  }
  fullPost[(nburn + 1):niter,]
}


################################################################################
### mcmcPS
################################################################################
mcmcPS_EKnown <- function(Estar, Istar, X, maxInf, S0, E0, I0, N, 
                   niter, nburn, inits,
                   dist, 
                   betaProVar, psParamsProVar,
                   betaPrior, rateEPrior, psParamsPrior,
                   WAIC, progress) {
  
  # useful things from data
  epiSize <- sum(Istar)
  tau <- NROW(Istar)
  popSizesInv <- N^(-1)
  
  # initialize path matrix through infectious period 
  infTime <- rep( which(Istar>0), Istar[ which(Istar>0)])
  lengthInfVec <- rep(round(maxInf)/2, length(infTime))
  fullI <- getFullIPS(tau, maxInf, infTime, lengthInfVec) 
  fullRstar <-getFullRStar(tau, maxInf, infTime, lengthInfVec)
  
  # initalize compartment vectors
  susVec <- getS(S0, Estar)
  expVec <- getE(E0, Estar, Istar)
  infVec <- rowSums(fullI)
  
  psParamNames <- names(inits$psParams)
  
  # initialize storage
  # variables in transmission probability (betas)
  transProbPost <- matrix(NA, nrow = niter, ncol = ncol(X))
  rateEPost <- rep(NA, niter)
  psParamPost <- matrix(NA, nrow = niter, ncol = length(psParamNames))
  colnames(psParamPost) <- psParamNames
  
  transProbPost[1,] <- inits$beta
  rateEPost[1] <- inits$rateE
  psParamPost[1,] <- inits$psParams
  
  # initialize covariance for transmission probability proposal
  S <- diag(betaProVar, ncol(transProbPost))
  
  # initialize storage of log-likelihood if WAIC is requested
  if (WAIC) {
    ll_store <- matrix(NA, nrow = length(Istar), ncol = niter)
    ll_store[,1] <- ll_PS(susVec=susVec, Estar=Estar, X=X, beta=transProbPost[1,], 
                          infVec=infVec, popSizesInv=popSizesInv, 
                          expVec=expVec, Istar=Istar, rateE=rateEPost[1],
                          maxInf=maxInf, dist=dist, psParams=psParamPost[1,], 
                          fullI=fullI, fullRstar=fullRstar)
  }
  
  # create progress bar
  if (progress) {
    pb <- txtProgressBar(min = 1, max = niter, style = 3)
  }
  
  for (i in 2:niter) {
    
    # transmission probability
    transProbUpdate <- MH_update_adapt(x0=transProbPost[i - 1, ], f=fcTransProb, 
                                       S=S, currentIter=i, nBurn=nburn,
                                       susVec=susVec, Estar=Estar, X=X,  
                                       infVec=infVec, popSizesInv=popSizesInv,
                                       betaPrior=betaPrior)
    
    transProbPost[i,] <- transProbUpdate$x
    S <- transProbUpdate$S
    
    # exponential rate parameter for E to I
    rateEPost[i] <- MH_update(rateEPost[i-1], fcRate, g_positive, g_s_positive, sd=0.05, 
                              compVec=expVec, transVec=Istar, ratePrior=rateEPrior)
    
    # path-specific parameters
    psParamPost[i,] <- MH_update(psParamPost[i-1,], fcPS, g_positive, g_s_positive,
                                 sd=psParamsProVar, 
                                 psParamNames=psParamNames, fullI=fullI, 
                                 fullRstar=fullRstar, dist=dist,
                                 maxInf=maxInf, psParamsPrior=psParamsPrior)
    
    # Rstar
    RstarProp <- update_Rstar_PS(lengthInf=lengthInfVec, epiSize=epiSize,
                                 maxInf=maxInf, infTime=infTime, fullI=fullI,
                                 fullRstar=fullRstar, infVec=infVec, 
                                 susVec=susVec, Estar=Estar, X=X,
                                 beta=transProbPost[i,], popSizesInv=popSizesInv,
                                 dist=dist, psParams=psParamPost[i,], percent=0.1)
    
    
    fullI <- RstarProp$fullI
    fullRstar <- RstarProp$fullRstar
    infVec <- RstarProp$infVec
    lengthInfVec <- RstarProp$lengthInf
    
    # update log-likelihood
    if (WAIC) {
      ll_store[,i] <- ll_PS(susVec=susVec, Estar=Estar, X=X, beta=transProbPost[i,], 
                            infVec=infVec, popSizesInv=popSizesInv, 
                            expVec=expVec, Istar=Istar, rateE=rateEPost[i],
                            maxInf=maxInf, dist=dist, psParams=psParamPost[i,], 
                            fullI=fullI, fullRstar=fullRstar)
    }
    
    # update progress bar
    if (progress) {
      setTxtProgressBar(pb, i)
    }
    
  }
  
  # export posteriors in one data frame
  fullPost <- cbind.data.frame(transProbPost, rateEPost, psParamPost)
  colnames(fullPost) <- c(paste0('beta', 1:ncol(transProbPost)),
                          'rateE', psParamNames)
  
  if (WAIC) {
    return(list(fullPost = fullPost[(nburn + 1):niter,], 
                WAIC = getWAIC(ll_store[,(nburn + 1):niter])))
  }
  fullPost[(nburn + 1):niter,]
}




################################################################################
### mcmcIDD
################################################################################
mcmcIDD_EKnown <- function(Estar, Istar, X, maxInf, S0, E0, I0, N, 
                    niter, nburn, inits,
                    iddFun, iddParamNames,
                    betaProVar, iddParamsProVar,
                    betaPrior, rateEPrior, iddParamsPrior,
                    WAIC, progress) {
  
  # useful things from data
  epiSize <- sum(Istar)
  tau <- NROW(Istar)
  popSizesInv <- N^(-1)
  
  # initialize matrix of infectious duration over epidemic time
  infTime <- rep( which(Istar>0), Istar[ which(Istar>0)])
  fullI <- getFullI(tau, maxInf, infTime) 
  
  # initalize compartment vectors
  susVec <- getS(S0, Estar)
  expVec <- getE(E0, Estar, Istar)
  infVec <- rowSums(fullI)
  
  # beta and idd param names
  betaNames <- paste0('beta', 1:ncol(X))
  
  # initialize storage
  # variables in transmission probability
  transProbPost <- matrix(NA, nrow = niter, 
                          ncol = length(betaNames) + length(iddParamNames))
  colnames(transProbPost) <- c(betaNames, iddParamNames)
  rateEPost <- rep(NA, niter)
  
  transProbPost[1,] <- unlist(c(inits$beta, inits$iddParams))
  rateEPost[1] <- inits$rateE
  
  # initialize covariance for transmission probability proposal
  S <- diag(c(betaProVar, iddParamsProVar), ncol(transProbPost))
  
  # initialize storage of log-likelihood if WAIC is requested
  if (WAIC) {
    ll_store <- matrix(NA, nrow = length(Istar), ncol = niter)
    ll_store[,1] <- ll_IDD(susVec=susVec, Estar=Estar, X=X, 
                           beta=transProbPost[1, betaNames],
                           popSizesInv=popSizesInv, iddFun=iddFun, 
                           iddParams=transProbPost[1, iddParamNames], 
                           fullI=fullI, maxInf=maxInf,
                           expVec=expVec, Istar=Istar, rateE=rateEPost[1]) 
  }
  
  # create progress bar
  if (progress) {
    pb <- txtProgressBar(min = 1, max = niter, style = 3)
  }
  
  for (i in 2:niter) {
    
    # transmission probability
    transProbUpdate <- MH_update_adapt(transProbPost[i - 1, ], fcTransProbIDD, 
                                       S, i, nburn,
                                       susVec=susVec, Estar=Estar, X=X,  
                                       popSizesInv=popSizesInv, iddFun=iddFun,
                                       betaNames=betaNames, iddParamNames=iddParamNames, 
                                       fullI=fullI, maxInf=maxInf, 
                                       betaPrior=betaPrior, iddParamsPrior=iddParamsPrior)
    
    transProbPost[i,] <- transProbUpdate$x
    S <- transProbUpdate$S
    
    # exponential rate parameter for E to I
    rateEPost[i] <- MH_update(rateEPost[i-1], fcRate, g_positive, g_s_positive, sd=0.05, 
                              compVec=expVec, transVec=Istar, ratePrior=rateEPrior)
    
   
    
    # update log-likelihood
    if (WAIC) {
      ll_store[,i] <- ll_IDD(susVec=susVec, Estar=Estar, X=X, 
                             beta=transProbPost[i, betaNames],
                             popSizesInv=popSizesInv, iddFun=iddFun, 
                             iddParams=transProbPost[i, iddParamNames], 
                             fullI=fullI, maxInf=maxInf,
                             expVec=expVec, Istar=Istar, rateE=rateEPost[i]) 
    }
    
    
    # update progress bar
    if (progress) {
      setTxtProgressBar(pb, i)
    }
    
  }
  
  # export posteriors in one data frame
  fullPost <- cbind.data.frame(transProbPost, rateEPost)
  colnames(fullPost) <- c(betaNames, iddParamNames, 'rateE')
  
  if (WAIC) {
    return(list(fullPost = fullPost[(nburn + 1):niter,], 
                WAIC = getWAIC(ll_store[,(nburn + 1):niter])))
  }
  fullPost[(nburn + 1):niter,]
}

