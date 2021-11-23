################################################################################
# Miscellaneous helper functions
################################################################################

# vectorized truncated distribution function for path-specific probability vectors
ptruncVec <- function(x, spec, a = -Inf, b = Inf, ...){
  tt <- x
  aa <- a
  bb <- rep(b, length(x))
  G <- get(paste("p", spec, sep = ""), mode = "function")
  tt <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
  tt <- tt - G(aa, ...)
  tt <- tt/(G(bb, ...) - G(aa, ...))
  return(tt)
}

# used in mcmcIDD to initialize the matrix off infectious durations over epidemic 
getFullI <- function(nTime, maxInf, infTime) {
  
  emptyMat <- matrix(0,nrow = nTime, ncol = maxInf)
  
  for (k in 1:length(infTime)) {
    infTime_i <- infTime[k]
    length_i <- maxInf
    
    rowIdx <- infTime_i + 1:length_i
    
    rowIdx <- rowIdx[rowIdx>0]
    colIdx <- 1:length_i
    colIdx <- colIdx[1:length(rowIdx)]
    
    emptyMat[rowIdx,colIdx] <- emptyMat[rowIdx,colIdx] + diag(length(rowIdx))
  }
  
  emptyMat
}

# used in mcmcPS to initialize the path-specific matrix
getFullIPS <- function(nTime, maxInf, infTime, lengthInf) {
  
  emptyMat <- matrix(0,nrow = nTime, ncol = maxInf)
  
  for (k in 1:length(infTime)) {
    infTime_i <- infTime[k]
    length_i <- lengthInf[k]
    
    rowIdx <- infTime_i + 1:length_i
    
    rowIdx <- rowIdx[rowIdx>0]
    colIdx <- 1:length_i
    colIdx <- colIdx[1:length(rowIdx)]
    
    emptyMat[rowIdx,colIdx] <- emptyMat[rowIdx,colIdx] + diag(length(rowIdx))
  }
  
  emptyMat
}

# used in mcmcPS to initialize the path-specific Rstar
getFullRStar <- function(nTime, maxInf, infTime, lengthInf) {
  
  emptyMat <- matrix(0,nrow = nTime, ncol = maxInf)
  for (k in 1:length(infTime)) {
    infTime_i <- infTime[k]
    length_i <- lengthInf[k]
    
    rowIdx <- infTime_i+length_i
    
    
    emptyMat[rowIdx,length_i] <- emptyMat[rowIdx,length_i] + 1
    
  }
  emptyMat
}


# update full I for one person (called in update_Rstar_PS)
updateFullIOne <- function(fullI, updateInfTime,
                           updateLengthInf, updateLengthInfNew) {
  fullNew <- fullI
  
  rowIdx <- updateInfTime + 1:updateLengthInf
  rowIdxNew <- updateInfTime + 1:updateLengthInfNew
  rowIdxNew <- rowIdxNew[rowIdxNew>0]
  
  colIdx <- 1:updateLengthInf
  colIdxNew <- 1:updateLengthInfNew
  colIdxNew <- colIdxNew[1:length(rowIdxNew)]
  
  # remove previous path
  if (updateLengthInf == 1) {
    fullNew[rowIdx,colIdx] <- fullNew[rowIdx,colIdx]- 1
  } else {
    diag(fullNew[rowIdx,colIdx]) <- diag(fullNew[rowIdx,colIdx]) - 1
  }
  
  # add new path
  if (updateLengthInfNew == 1) {
    fullNew[rowIdxNew,colIdxNew] <- fullNew[rowIdxNew,colIdxNew] + 1
  } else {
    diag(fullNew[rowIdxNew,colIdxNew]) <- diag(fullNew[rowIdxNew,colIdxNew]) + 1
  }
  
  fullNew
}

# update full Rstar for one person (called in update_Rstar_PS)
updateFullRstarOne <- function(fullRstar, updateInfTime,
                               updateLengthInf, updateLengthInfNew) {
  fullNew <- fullRstar
  
  rowIdx <- updateInfTime + updateLengthInf
  rowIdxNew <- updateInfTime + updateLengthInfNew
  
  fullNew[rowIdx,updateLengthInf] <- fullNew[rowIdx,updateLengthInf] - 1
  fullNew[rowIdxNew,updateLengthInfNew] <- fullNew[rowIdxNew,updateLengthInfNew] + 1
  fullNew
}



#' Get path-specific vector of removal probabilities over duration of infection
#'
#' @param maxInf Integer. The maximum length of time an individual can 
#' spend in the infectious compartment
#' @param dist Character. 
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' @export
psProbVec <- function(maxInf, dist, psParams) {
  timeInf <- 1:maxInf
  nextTime <- timeInf + 1
  
  if (dist == 'exp') {
    transProb <- ptruncVec(nextTime, spec = dist, 
                           a = timeInf, b = Inf, 
                           rate = psParams$rate)
  } else if (dist == 'gamma') {
    transProb <- ptruncVec(nextTime, spec = dist, 
                           a = timeInf, b = Inf, 
                           shape = psParams$shape, rate = psParams$rate)
  } else if (dist == 'weibull') {
    transProb <- ptruncVec(nextTime, spec = dist, 
                           a = timeInf, b = Inf, 
                           shape = psParams$shape, scale = psParams$scale)
  }
  
  
  transProb[maxInf] <- 1
  transProb[is.na(transProb)] <- 1
  
  transProb
}


psParamsCheck <- function(psParams, dist) {
  if (dist == 'exp' & !all(names(psParams) %in% c('rate'))) {
    stop("if dist = 'exp', psParams must contain rate parameter")
    
    psParams <- psParams
    
  } else if (dist == 'gamma') {
    if (!all(names(psParams) %in% c('shape', 'rate'))) {
      if (!all(names(psParams) %in% c('shape', 'scale'))) {
        stop("if dist = 'gamma', psParams must contain shape and rate/scale parameters")
      }
      psParams <- list(shape = psParams$shape,
                    rate = 1 / psParams$scale)
      names(psParams) <- c('shape', 'rate')
    }
    # make sure order is correct
    psParams <- list(shape = psParams$shape,
                  rate = psParams$rate)
    names(psParams) <- c('shape', 'rate')
    
  } else if (dist == 'weibull') {
    if (!all(names(psParams) %in% c('shape', 'scale'))) {
      stop("if dist = 'weibull', psParams must contain shape and scale parameters")
    }
    # make sure order is correct
    psParams <- list(shape = psParams$shape,
                  scale = psParams$scale)
    names(psParams) <- c('shape', 'scale')
  }
  psParams
}


getWAIC <- function(ll) {
  
  pWAIC2 <- sum(apply(ll, 1, var))
  
  lppd <- sum(log(rowMeans(exp(ll))))
  
  -2 * (lppd - pWAIC2)
  
}


# fix args
fixFunArgs <- function(x, env = parent.frame()) {
  #x <- substitute(x, env = env) # substitute happens in bootstrap
  fname <- deparse1(x[[1]])
  arggs <- as.list(x)[-1]
  w <- get(fname)
  for (nn in names(arggs)) {
    formals(w)[[nn]] <- arggs[[nn]]
  }
  w
}
