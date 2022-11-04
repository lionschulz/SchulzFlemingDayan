# This file contains the function which computes the summary statistics for the
# Rmarkdowns. It uses numerical integration over predefined grids to achieve this.

# The code is organized as follows
# 1 - Basics:       defines the grid
# 2 - Computations: defines the model
# 3 - Statistics:   gets the statistics
# 4 - Output:       defines what is returned

# The statistics are computed by integrating over the probability space between
# specific intervals. For example, for the average seeking, we find the points
# where the two Q-value curves intersect and then take the probability mass 
# between the intersections as the probability of search.

# Apart from the model parameters (zetaI,sF,rS), the functions also contains two
# further arguments
# - aF:   If TRUE, the function also computes the final accuracy statistics.
#         This can be switched off if we want the code to run faster and don't
#         these statistics
# - resolution: Defines the resolution of the grid used for numerical integration.


model.Pd <-  function(zetaI,# zeta_I  (joint standard dev. of actor and rater signal)
                      sF,   # sigma_F (standard deviation of final signal)
                      rS,   # r_S     (cost of information seeking)
                      aF = T,
                      resolution = 2000){
  
  status_message <- paste("Postdecisional model with: zeta_I =", round(zetaI, 3),
                          ", cost = ", rS,
                          ", sigma_F = ", sF,
                          ", aF setting = ", aF)
  
  print(status_message)
  
  # Basic Set-Up ---- ##########################################################
  d <- 1
  space_limits <- qnorm(10e-7,0,max(zetaI,sF)) # limits space for 
  
  # this vector forms the basis for the grids representing the variable space 
  # further down below. It is used for all random variables (Z_I and X_F)
  x <- seq(from = space_limits, to = -space_limits, length.out = resolution) 
  
  
  ## Computations ---- #########################################################
  
  # Q_S(0) ---- (eq. 11)
  posterior_ZI <- sigmoid(2*x/zetaI^2) # equiv. to computing the confidence (eq. 6)
  Q0 <- pmax(posterior_ZI, 1 - posterior_ZI)
  
  
  ## Q(1) ---- 
  zetaF2 <- joinsd(zetaI,sF)^2  # getting the combined variance, (eq. 12)
  ZF <- outer(x/zetaI^2,x/sF^2, FUN = "+")*zetaF2 # creating a grid of ZF
  posterior_ZF <- sigmoid(2*ZF/zetaF2) # getting the posterior for it
  value_ZF <- 0.5 + abs(posterior_ZF -.5) # max of the two posteriors (eq. 13)
  

  # P(XF|ZI) ----
  # Setting up the likelihoods for P(XF|ZI) (see eq. 22)
  p_XF <- dnorm(x,mean = d,sd = sF)
  p_XF <- normalize(p_XF)
  p_XF_minus <- dnorm(x,mean = -d,sd = sF)
  p_XF_minus <- normalize(p_XF_minus)
  p_XFZI <- outer(posterior_ZI,p_XF) + outer(1-posterior_ZI,p_XF_minus)
  
  Q1 <- rowSums(p_XFZI*value_ZF) + rS # (see eq. 14)
  
  
  ## Statistics ---- ###########################################################
  
  ## Average seeking ----
  
  # Getting the difference between the two Q-values
  Q_difference <- Q1 - Q0 
  Q_difference <- sign(Q_difference)
  # Identify the intersection (only necessary for half the space b/c of symmetry)
  diffcheck <- diff(Q_difference[1:(resolution/2)]) # returns 2 at intersection
  
  # compute the average seeking
  if (rS == 0 ){ # seeking is free 
    # this avoids numerical errors that can appear when seeking is free 
    # -> seeking is always valuable, then, but the Q-value curves might be very close
    seek_zone_cutoff <- -100
    av_search <- 1 # resolution_grid
    
  } else if (min(Q_difference) >= 0) { 
    # seeking is always as or more valuable than not seeking
    seek_zone_cutoff <- -100
    av_search <- 1
    
  } else if (max(Q_difference) == -1){ 
    # seeking is always less valuable than not seeking
    seek_zone_cutoff <- 0
    av_search <- 0
    
  } else {
    # regular seeking case: the Q-values intersect somewhere
    index_seek <- which(diffcheck != 0)
    # we get the random variable value where the two intersect
    seek_zone_cutoff <- mean(x[index_seek:(index_seek+1)])
    # take the mass between the two points
    av_search <- pnorm(-1*seek_zone_cutoff,d,zetaI) - pnorm(seek_zone_cutoff,d,zetaI)
  }
  
  
  ## Average final correctness ----
  
  if (aF  ==  T) { # If we want to compute it:
    
    # Index matrix where sought and where correct (sign(ZF = d))
    index_sought_correct <- ZF > 0 & Q_difference == 1
    index_sought_incorrect <- ZF < 0 & Q_difference == 1

    # Get densities for ZF (two dimensional)
    p_XI <- dnorm(x,d,zetaI)
    p_ZF <- outer(p_XI, dnorm(x,d,sF))
    
    # Accuracy after seeking: sum the densities where sought
    av_accuracy_s1 <- sum(index_sought_correct*p_ZF)
    incorrect_s1 <- sum(index_sought_incorrect*p_ZF)
    av_accuracy_s1 <- av_accuracy_s1/(av_accuracy_s1+incorrect_s1)
    
    # Accuracy after not seeking: ratio of positve and negative densities
    av_accuracy_s0 <- 1 - pnorm(-seek_zone_cutoff,1,zetaI) # sum(seekZvec_nc*xdens)
    incorrect_s0 <- pnorm(seek_zone_cutoff,1,zetaI)
    av_accuracy_s0 <- av_accuracy_s0/(av_accuracy_s0+incorrect_s0)
    
  } else { # if we don't want to compute all the final accuracy averages
    av_accuracy_s1 <- NA
    av_accuracy_s0 <- NA
  }
  
  
  ## Confidence cut-off ----
  
  if (av_search != 1 & av_search != 0) {
    # When seeking != 100 % or 0 %
    # take the random variable where Q0 = Q1 and compute the confidence
    cI_cutoff <- sigmoid(-2*seek_zone_cutoff/zetaI^2) # 
  } else {
    cI_cutoff <- NA
  }
  
  
  ## Output ----- ##############################################################
  
  output <- list()
  output$av_search <- av_search
  output$av_accuracy_s1 <- av_accuracy_s1
  output$av_accuracy_s0 <- av_accuracy_s0
  output$cI_cutoff <- cI_cutoff
  return(output)
}  
