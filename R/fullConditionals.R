################################################################################
# compiled functions for full conditionals
################################################################################

### full conditional for parameters in the transmission probability
# for exponential and path specific, this is just beta's

fcTransProb <- function(x, susVec, Estar, X, infVec, popSizesInv, betaPrior) {
  
  f_SE(susVec=susVec, Estar=Estar, X=X, beta=x, 
       infVec=infVec, popSizesInv=popSizesInv) +
    betaPrior(x)
}

fcTransProbIDD <- function(x, susVec, Estar, X, popSizesInv, 
                           iddFun, betaNames, iddParamNames, fullI, maxInf,
                           betaPrior, iddParamsPrior) {
  
  # split x into iddParams and beta
  iddParams <- x[iddParamNames] 
  beta <- x[betaNames]
  
  f_SE_IDD(susVec=susVec, Estar=Estar, X=X, beta=beta, popSizesInv=popSizesInv, 
           iddFun=iddFun, iddParams=iddParams, fullI=fullI, maxInf=maxInf) + 
    betaPrior(beta) + 
    iddParamsPrior(iddParams)
}


### full conditional for rate parameter from exponential transitions
fcRate <- function(x, compVec, transVec, ratePrior) {
  f_exp(compVec=compVec, transVec=transVec, rate=x) + 
    ratePrior(x)
}

### full conditional for path-specific parameters (proposed jointly)
fcPS <- function(x, psParamNames, fullI, fullRstar, dist, maxInf, psParamsPrior) {
  names(x) <- psParamNames
  
  f_PS(fullMat=fullI, fullStar=fullRstar, dist=dist, 
       psParams=x, maxLength=maxInf) + 
    psParamsPrior(x)
}



### full conditional for Estar vector 
# exponential and PS
fcEstar <- function(susVec, Estar, X, beta, infVec, popSizesInv,
                    expVec, Istar, rateE) {
  
  f_SE(susVec=susVec, Estar=Estar, X=X, beta=beta, 
       infVec=infVec, popSizesInv=popSizesInv) +
    f_exp(compVec=expVec, transVec=Istar, rate=rateE)
  
}

fcEstar_IDD <- function(susVec, Estar, X, beta, popSizesInv,
                    iddFun, iddParams, fullI, maxInf,
                    expVec, Istar, rateE) {
  
  f_SE_IDD(susVec=susVec, Estar=Estar, X=X, beta=beta, popSizesInv=popSizesInv, 
           iddFun=iddFun, iddParams=iddParams, fullI=fullI, maxInf=maxInf) +
    f_exp(compVec=expVec, transVec=Istar, rate=rateE)
  
}

###  full conditional for Rstar
# exponential transition
fcRstar <- function(susVec, Estar, X, beta, infVec, popSizesInv,
                    Rstar, rateI) {
  
  # infVec changed
  f_SE(susVec=susVec, Estar=Estar, X=X, beta=beta, 
       infVec=infVec, popSizesInv=popSizesInv)  + 
    f_exp(compVec=infVec, transVec=Rstar, rate=rateI)
}


# path-specific transition
fcRstar_PS <- function(fullI, fullRstar, susVec, Estar, X, beta,
                    infVec, popSizesInv, dist, psParams, maxInf) {
  
  # infVec, fullI, fullRstar also changes
  f_SE(susVec=susVec, Estar=Estar, X=X, beta=beta, 
       infVec=infVec, popSizesInv=popSizesInv) + 
    f_PS(fullMat=fullI, fullStar=fullRstar, dist=dist, 
         psParams=psParams, maxLength=maxInf) 
  
}


