
################################################################################
# compiled functions to update transition compartments
# Exponential transition from E to I
# Exponential transition from I to R
# Path specific transition from I to R
################################################################################

### Exponential transition from E to I 
# Path-specific and exponential model (constant transmission)

update_Estar <- function(Estar, S0, E0, Istar, susVec, expVec, 
                         X, beta, infVec, popSizesInv, rateE, 
                         percent = 0.1) {
  
  # iterate the following procedure a number of times (10% of final size)
  n_iter <- floor(percent*sum(Estar))
  
  for (i in 1:n_iter) {
    # find candidate that satisfies all constraints
    it <- 0
    repeat {
      EstarNew <- Estar
      
      it <- it + 1
      
      # select a time point  with a non-zero value to subtract 1 from
      if (length(which(Estar>0))==1) {
        positiveidx <- sample(which(Estar>0), 1,prob = c(rep(0, which(Estar>0) - 1), 1))
      } else {
        positiveidx <- sample(which(Estar>0), 1)
      }
      EstarNew[positiveidx] <- Estar[positiveidx] -1
      
      # select any time point to add 1 to
      # cannot become exposed after the last Istar
      idx <- sample(1:max(which(Istar > 0)), 1)
      EstarNew[idx] <- EstarNew[idx] + 1
      
      # update S, E
      susVecNew <- getS(S0, EstarNew)
      expVecNew <- getE(E0, EstarNew, Istar)
      
      # check conditions (E(t) >= 0, E(t) + I(t) > 0, E(tau) + I(tau) = 0)
      condition1 <- sum(expVecNew<0)==0
      condition2 <- sum(head(expVecNew,-1)+head(infVec,-1) < 0)==0
      
      # if not true try again
      if ((condition1&condition2)==TRUE) {
        break
      } else if (it >= 1000) {
        break
      }
    }
    
    # accept EstarNew with probability min(1, f(EstarNew)/f(Estar)) 
    # on log scale this is probability min(0, log f(EstarNew) - log f(Estar))
    f1 <- fcEstar(susVec=susVecNew, Estar=EstarNew, X=X, beta=beta, 
                  infVec=infVec, popSizesInv=popSizesInv,
                  expVec=expVecNew, Istar=Istar, rateE=rateE)
    
    f0 <- fcEstar(susVec=susVec, Estar=Estar, X=X, beta=beta, 
                  infVec=infVec, popSizesInv=popSizesInv,
                  expVec=expVec, Istar=Istar, rateE=rateE)
    
    a <- min(0, f1 - f0)
    u <- log(runif(1, 0, 1))
    
    if (u < a){
      # accept
      Estar <- EstarNew
      susVec <- susVecNew
      expVec <- expVecNew
    }
    # else everything remains as is
    
  }
  
  return(list(Estar = Estar, susVec = susVec, expVec=expVec))
  
}


update_Estar_IDD <- function(Estar, S0, E0, Istar, susVec, expVec, infVec,
                             X, beta, popSizesInv, 
                             iddFun, iddParams, fullI, maxInf, rateE, 
                             percent = 0.1) {
  
  # iterate the following procedure a number of times (10% of final size)
  n_iter <- floor(percent*sum(Estar))
  
  for (i in 1:n_iter) {
    # find candidate that satisfies all constraints
    it <- 0
    repeat {
      EstarNew <- Estar
      
      it <- it + 1
      
      # select a time point  with a non-zero value to subtract 1 from
      if (length(which(Estar>0))==1) {
        positiveidx <- sample(which(Estar>0), 1,prob = c(rep(0, which(Estar>0) - 1), 1))
      } else {
        positiveidx <- sample(which(Estar>0), 1)
      }
      EstarNew[positiveidx] <- Estar[positiveidx] -1
      
      # select any time point to add 1 to
      # cannot become exposed after the last Istar
      idx <- sample(1:max(which(Istar > 0)), 1)
      EstarNew[idx] <- EstarNew[idx] + 1
      
      # update S, E
      susVecNew <- getS(S0, EstarNew)
      expVecNew <- getE(E0, EstarNew, Istar)
      
      # check conditions (E(t) >= 0, E(t) + I(t) > 0, E(tau) + I(tau) = 0)
      condition1 <- sum(expVecNew<0)==0
      condition2 <- sum(head(expVecNew,-1)+head(infVec,-1) < 0)==0
      
      # if not true try again
      if ((condition1&condition2)==TRUE) {
        break
      } else if (it >= 1000) {
        break
      }
    }
    
    # accept EstarNew with probability min(1, f(EstarNew)/f(Estar)) 
    # on log scale this is probability min(0, log f(EstarNew) - log f(Estar))
    f1 <- fcEstar_IDD(susVec=susVecNew, Estar=EstarNew, X=X, beta=beta, 
                      popSizesInv=popSizesInv, iddFun=iddFun, 
                      iddParams=iddParams, fullI=fullI, maxInf=maxInf,
                      expVec=expVecNew, Istar=Istar, rateE=rateE)
    
    f0 <- fcEstar_IDD(susVec=susVec, Estar=Estar, X=X, beta=beta, 
                      popSizesInv=popSizesInv, iddFun=iddFun, 
                      iddParams=iddParams, fullI=fullI, maxInf=maxInf,
                      expVec=expVec, Istar=Istar, rateE=rateE)
    
    a <- min(0, f1 - f0)
    u <- log(runif(1, 0, 1))
    
    if (u < a){
      # accept
      Estar <- EstarNew
      susVec <- susVecNew
      expVec <- expVecNew
    }
    # else everything remains as is
    
  }
  
  return(list(Estar = Estar, susVec = susVec, expVec=expVec))
  
}


### Exponential transition from I to R 
update_Rstar_exp <- function(Rstar,  maxInf, epiSize, 
                             susVec, Estar, X, beta,
                             infVec, popSizesInv, rateI,
                             expVec, I0, Istar, percent) {
  
  
  fOld <- fcRstar(susVec=susVec, Estar=Estar, X=X, beta=beta,
                  infVec=infVec, popSizesInv=popSizesInv, 
                  Rstar=Rstar, rateI=rateI)
  
  lastRemTime <- max(which(Istar > 0)) + maxInf
  firstInfTime <- min(which(Istar > 0))  + 1
  
  #number updates to do
  nUpdate <- min(floor(percent * epiSize), 5000)
  
  updateCount <- 0
  
  repeat{
    
    # select a time point  with a non-zero value to subtract 1 from
    toUpdateIdx <- sample(which(Rstar > 0), 1)
    
    RstarNew <- Rstar
    
    RstarNew[toUpdateIdx] <- Rstar[toUpdateIdx] - 1
    
    toAddIdx <- sample(firstInfTime:lastRemTime, 1)
    RstarNew[toAddIdx] <- RstarNew[toAddIdx] + 1
    
    # update I, R
    infVecNew <- getI(I0, Istar, RstarNew)
    
    # check conditions (E(t) + I(t) > 0, E(tau) + I(tau) = 0, I(t) >=0)
    lastInfTime <- max(which(infVecNew > 0))
    
    condition1 <- sum(expVec[1:lastInfTime] + infVecNew[1:lastInfTime] < 0)==0
    condition2 <- tail(expVec,1) + tail(infVecNew,1) == 0
    condition3 <- sum(infVecNew<0)==0
    
    if (condition1 & condition2 & condition3) {
      
      # calculate full conditional for new 
      # old has already been calculated after the first iteration
      fNew <- fcRstar(susVec=susVec, Estar=Estar, X=X, beta=beta,
                      infVec=infVecNew, popSizesInv=popSizesInv,
                      Rstar=RstarNew, rateI=rateI)
      
      a <- fNew - fOld
      a <- min(0, a)
      
      u <- log(runif(1,0,1))
      
      
      if (u < a) {
        
        Rstar <- RstarNew
        infVec <- infVecNew
        fOld <- fNew
        updateCount <- updateCount + 1
        
      } 
      
    }
    
    if (updateCount == nUpdate) break
  }
  
  return(list(Rstar=Rstar, infVec=infVec))
  
}



# Path specific transition from I to R
update_Rstar_PS <- function(lengthInf, epiSize, maxInf,
                            infTime, fullI, fullRstar, infVec,
                            susVec, Estar, X, beta,
                            popSizesInv, dist, psParams, percent) {
  
  
  fOld <- fcRstar_PS(fullI=fullI, fullRstar=fullRstar, 
                     susVec=susVec, Estar=Estar, X=X, beta=beta,
                     infVec=infVec, popSizesInv=popSizesInv,
                     dist=dist, psParams=psParams, maxInf=maxInf)
  
  #number updates to do
  nUpdate <- min(floor(percent * epiSize), 5000)
  
  updateCount <- 0
  
  repeat{
    
    # randomly select an individual
    toUpdateIdx <- sample(1:length(lengthInf), 1)
    
    lengthInfNew <- lengthInf
    lengthInfNew[toUpdateIdx] <- sample(1:maxInf, 1)
    
    # check that proposed time was during the epidemic
    condition1 <- infTime[toUpdateIdx] + lengthInfNew[toUpdateIdx] <= NROW(susVec)
    
    
    if (condition1 & (lengthInfNew[toUpdateIdx] != lengthInf[toUpdateIdx])) {
      
      fullINew <- updateFullIOne(fullI, infTime[toUpdateIdx],
                                 lengthInf[toUpdateIdx], lengthInfNew[toUpdateIdx])
      fullRstarNew <- updateFullRstarOne(fullRstar,  infTime[toUpdateIdx],
                                         lengthInf[toUpdateIdx], lengthInfNew[toUpdateIdx])
      
      infVecNew <-  rowSums(fullINew)
      
      # calculate full conditional for new 
      # old has already been calculated after the first iteration
      fNew <- fcRstar_PS(fullI=fullINew, fullRstar=fullRstarNew, 
                         susVec=susVec, Estar=Estar, X=X, beta=beta,
                         infVec=infVecNew, popSizesInv=popSizesInv,
                         dist=dist, psParams=psParams, maxInf=maxInf)
      
      a <- fNew - fOld
      a <- min(0, a)
      
      u <- log(runif(1,0,1))
      
      
      if (u < a) {
        
        lengthInf <- lengthInfNew
        fullI <- fullINew
        fullRstar <- fullRstarNew
        infVec <- infVecNew
        fOld <- fNew
        
        updateCount <- updateCount + 1
        
      } 
      
    }
    
    if (updateCount == nUpdate) break
  }
  
  return(list(fullI = fullI, lengthInf = lengthInf, 
              fullRstar=fullRstar, infVec=infVec))
  
}








