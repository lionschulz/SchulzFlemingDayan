# Functions for the second-order model

# Function that sets up the covariance function for the bivariate normal dens
computeCov = function(sd1, sd2, rho){
  CapitalSigma = matrix(c(sd1^2, rho*sd1*sd2, rho*sd1*sd2, sd2^2), ncol = 2)
}


# Second-order confidence
# Adapted from Matlab Code from Fleming & Daw, 2017
# https://github.com/smfleming/Self-evaluation-paper/blob/master/computeMetaConf.m
# Works on vectors of YI
# See appendix of paper for details on this.
cI.2o <- function(YI,sI,tI,rho,a,tol = 10e-100){

  # Initialization ----
  dhat <- c(-1, 1)        # the two actual d locations
  
  # Conditional distributions ----
  # mean XI|YI
  mu_cond <- matrix(NA,nrow = length(YI),ncol = length(dhat))    
  mu_cond[,1] <- dhat[1] + (sI/tI)*rho*(YI-dhat[1])
  mu_cond[,2] <- dhat[2] + (sI/tI)*rho*(YI-dhat[2])
  # standard deviation XI|YI
  sI_cond = sqrt((1-rho^2)*(sI^2))
  
  # p(XI|a)
  pXIa <- matrix(NA,nrow = length(YI),ncol = length(dhat))
  pXIa[a == -1,] <- pnorm(0, mu_cond[a == -1,], sI_cond)
  pXIa[a == 1,] <- pnorm(0, -mu_cond[a == 1,], sI_cond)
  
  # Likelihoods p(d|YI) ----
  liks <- matrix(NA,nrow = length(YI),ncol = length(dhat))    
  liks[,1] <- dnorm(YI, dhat[1], tI)
  liks[,2] <- dnorm(YI, dhat[2], tI)
  liks <- liks/rowSums(liks)
  
  # Checking tolerance ----
  liks[liks < tol] <- tol
  pXIa[pXIa < tol] <- tol
  
  # Normalization ----
  pfull <- liks*pXIa
  pfull <- pfull/apply(pfull, 1, FUN = sum)
  
  # Full output ---- 
  cI <- 1:length(YI)
  cI[a == -1] <- pfull[a == -1,1]
  cI[a == 1] <- pfull[a == 1,2]
  
  return(cI)
}



# Cue combination with correlations ----
# See appendix for details.
cue_combination.2o <- function(sI,tI,rho){
  
  output <- list()
  
  # Computations ----
  rXI <- 1/sI^2 # precIsion
  rYI <- 1/tI^2
  rcXI <- rXI - rho*sqrt(rXI*rYI) # the reliabilities including the correlations
  rcYI <- rYI - rho*sqrt(rXI*rYI)
  rcXYI <- rcXI + rcYI # the two together
  rZI <- (rXI + rYI - 2*rho*sqrt(rXI*rYI))/(1-rho^2)

  # Output ----
  output$wXI <- rcXI/rcXYI
  output$wYI <- 1 - output$wXI
  output$zetaI <- sqrt(1/rZI)
  return(output)
}

