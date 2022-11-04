# This file contains the function which computes the summary statistics for the
# Rmarkdowns. It uses numerical integration over predefined grids to achieve this.

# Note that this is the 'subjective' second-order model which is able to produce 
# both over and underconfidence. See 'model_2o.R' for the vanilla model.
# Also note that the subjective model is adapted not to do any final accuracy 
# calculations

# The code is organized as follows
# 1 - Basics:       defines the grid
# 2 - Computations: defines the model
# 3 - Statistics:   gets the statistics
# 4 - Output:       defines what is returned

# The statistics are computed by integrating over the probability space between
# specific intervals. For example, for the average seeking, we find the points
# where the two Q-value curves intersect and then take the probability mass 
# between the intersections as the probability of search.

# Apart from the model parameters (sI,sF,rS), the functions also contains two
# further arguments
# - output: what statistics are being put out
# - aF:     If TRUE, the function also computes the final accuracy statistics.
#           This can be switched off if we want the code to run faster and don't
#           these statistics
# - resolution: Defines the resolution of the grid used for numerical integration.
#           Note how this will be lower than the postdecisional model because
#           Integration in three-dimensional arrays is necessary, which places
#           increased demands on RAM
# - Qdifftol: Tolerance to avoid numerical problems when the two Q-values intersect
#           but never below a given threshold (defined by Qdifftol) which can 
#           give rise to problems in the computation of the statistics
# - denstol: tolerance for the lowest possible likelihood in the multivariate
#           normal distributions.

# We note that this model often only operates on one half of the X_I space
# which is possible due to the symmetry of the likelihoods and the confidence
# (which is only a function of Y_I and a_I)

model.2o_subjective <- function(
  sI_objective,   # sigma_I objective standard deviation of actor signal
  sI_subjective,  # subjective standard deviation of actor signal used by rater
  tI,   # tau_I      (standard deviation of rater signal)
  rho,  # rho_I      (correlation between actor/rater signal)
  sF,   # sigma_F    (standard deviation of final signal)
  rS,   # r_S        (cost of seeking)
  output = "av_search",
  aF = F, # compute final accuracy -> note that the subjective model has not been adjusted for this
  resolution = 200,
  Qdifftol = .001,
  denstol = 10e-300,
  subjective_confidence = FALSE) {
  
  # status print
  status_message <- paste("2nd-order model with: sigma_I_obj =", round(sI_objective, 3),
                          ", sigma_I_subj =", round(sI_subjective, 3),
                          ", tau_I =", round(tI, 3),
                          ", and rho = ", rho,
                          ", cost = ", rS,
                          ", sigma_F = ", sF,
                          ", aF setting = ", aF)
  
  print(status_message)
  
  # Basics ---- ################################################################
  
  # Check whether perfect calibration used
  sI_subjective <- ifelse(subjective_confidence, sI_subjective, sI_objective)
  
  d <- 1
  # compute the covariance matrix (see eq. 8)
  SI_objective <- computeCov(sI_objective, tI, rho) 
  SI_subjective <- computeCov(sI_subjective, tI, rho) 
  
  # checking the max sds to define the range of the grid
  maxsd <- max(sI_objective, tI, sF)
  mass_to_catch <- .9999 # defines the mass that we want to catch
  space_limits <- qnorm(mass_to_catch, 1, maxsd)
  
  # Setting up the grid using a beta distribution
  # This afford a higher resolution in the center while balancing numerical issues
  beta_parameter <- 2 # parameter defining the stretch of the beta distribution
  probability_space <-  seq(0,1,
                            length.out = resolution) # probability space to cover
  x <- (qbeta(probability_space, beta_parameter, beta_parameter)-.5)*
    (2*space_limits) # setting up the space of random variables
  
  # Getting the grid to work with the diffs
  xlong <- c(x[1]- 10e-6, x)
  diffx <- diff(xlong)
  diffxm <- outer(diffx,diffx)
  
  # Set up matrix with the densities of XI and YI
  # this is used for the predictions the agent makes and therefore uses 
  # the subjective confidence.
  xy <- expand.grid(x,x)
  p_XIYI_matrix <- matrix(dmvnorm(xy,
                                  mean = c(d,d),
                                  sigma = SI_subjective),
                          nrow = length(x))
  
  # Avoid underflow through the tolerance
  min_dens <- min(p_XIYI_matrix[p_XIYI_matrix > denstol])
  p_XIYI_matrix[p_XIYI_matrix < denstol] <- min_dens
  
  
  ## Computations ---- #########################################################
  
  # First-confidence ###########################################################
  
  # Computing the confidence (see eqs. 9 and 23 - 26) 
  cI_2o_ap <- cI.2o(x, sI_subjective, tI, rho, 1)     # confidence for positive values of xi
  cI_2o_am <- cI.2o(x, sI_subjective, tI, rho, -1)    # ... and for negative values
  
  # Getting this into matrices
  cI_mat_p <- matrix(cI_2o_ap, 
                     nrow = resolution/2, 
                     ncol = length(cI_2o_ap),
                     byrow = TRUE)
  cI_mat_m <- matrix(cI_2o_am, 
                     nrow = resolution/2,
                     ncol = length(cI_2o_am),
                     byrow = TRUE)
  
  
  # Seeking decision ###########################################################
  
  ## Linear weighting of the two cues taking into account the covariance ##
  # (see eqs. 31 - 36)
  cue_combination_list <- cue_combination.2o(sI_subjective, tI, rho) 
  zetaI <- cue_combination_list$zetaI
  wXI <- cue_combination_list$wXI
  wYI <- cue_combination_list$wYI
  
  
  ### Q-value for not seeking ## ----
  
  # (see eqs. 16, 17, appendix)
  
  ## Posterior based on ZI, and the associated value (eq. 37)
  ZI_matrix <- outer(x*wXI,x*wYI,FUN = "+")
  posterior_ZI <- sigmoid(2*ZI_matrix/zetaI^2) 
  posterior_ZIm <- 1 - posterior_ZI # the opposite of that posterior (for d = -1)
  # Value of ZI: 
  value_ZI <- pmax(posterior_ZI, posterior_ZIm) 
  
  
  ## Likelihood of XI given a YI and aI (eq. 40)
  # First, likelihood matrix of potential other source (d = -1)
  p_XIYI_matrix_minus <- matrix(dmvnorm(xy, 
                                        mean = c(-d,-d),
                                        sigma = SI_subjective),
                                nrow = length(x))
  # Then, we get half of the matrix of the positive source likelihood matrix
  p_XIYI_matrix_p <- p_XIYI_matrix[(dim(p_XIYI_matrix)[1]/2+1):
                                     dim(p_XIYI_matrix)[1],]
  # ... and half of the matrix of the negative source likelihood matrix
  p_XIYI_matrix_pm <- p_XIYI_matrix_minus[(dim(p_XIYI_matrix_minus)[1]/2+1): 
                                            dim(p_XIYI_matrix_minus)[1],]
  
  # We weight the two likelihood matrices set-up above 
  # according to the confidence
  p_XI_cond_YI <- p_XIYI_matrix_p * cI_mat_p + p_XIYI_matrix_pm * (1 - cI_mat_p) 
  p_XI_cond_YI_norm <- t(t(p_XI_cond_YI) / colSums(p_XI_cond_YI)) # making each row  sum to one
  
  # Getting half of the value
  value_ZI_ap <- value_ZI[(dim(value_ZI)[1]/2+1):dim(value_ZI)[1],]
  
  # Computing the value of not-seeking (see eq 16, 17)
  Q0 <- colSums(value_ZI_ap*p_XI_cond_YI_norm) # value for not seeking
  # this is now as a function of ZI, but only for one 'side' of XI
  # the other side of XI however has the same properties
  
  
  ### Q-value for seeking ## --- 
  
  # (see eq. 18-20, and appendix)
  
  # compute the overall standard deviation (no correlation, so normal function)
  zetaF <- joinsd(zetaI,sF) # (see eq. 12)
  
  ## Computing the value of ZF ##
  
  # We first do this using the likelihoods b/c this is faster than using sigmoid
  # However, it can produce 0 for both p_ZF_plus and p_ZF_minus which causes errors
  # In these cases value_ZF_cube ends up containing NAs - if this is the case,
  # we revert back to the sigmoid() method
  p_XF <- dnorm(x,d,sF)
  p_XF_minus <- dnorm(x,-d,sF)
  # Likelihoods of of the two possible sources in 3d space (this is all in aI = 1 space)
  p_ZF_plus <- outer(p_XIYI_matrix_p,p_XF) 
  p_ZF_minus <- outer(p_XIYI_matrix_pm,p_XF_minus) 
  posterior_ZF <- p_ZF_plus/(p_ZF_plus + p_ZF_minus)
  
  p_ZF_minus <- 1 - posterior_ZF
  value_ZF_cube <- pmax(posterior_ZF,p_ZF_minus) # 3d cube with values for ZF (eq. 41)
  
  # use sigmoid() method instead (as discussed above)
  if (sum(is.na(value_ZF_cube)) > 0) {
    
    ## Cue combination
    # in 3d, where we again just just consider half of the space
    ZI_matrixhalf <- ZI_matrix[(dim(ZI_matrix)[1]/2+1):dim(ZI_matrix)[1],]
    ZFcube <- outer(ZI_matrixhalf/zetaI^2, x/sF^2, FUN = "+")*zetaF^2
    # dim(ZFcube)
    
    # Getting the posterior for ZF
    posterior_ZF <- sigmoid(2*ZFcube/zetaF^2)
    # dim(posterior_ZF)
    posterior_ZFm <- 1 - posterior_ZF
    posterior_ZFm[is.na(posterior_ZFm)] <- 0 # avoids underflow
    value_ZF_cube <- pmax(posterior_ZF,posterior_ZFm)
  }
  
  # ... we now have the value but still need to predict XF
  
  ## Distribution of XF and XI as a function of YI and aI (eq. 44)
  # Set up distribution of XF as a function of YI (x: YI, y: XF, z: reps)
  p_XF_cond_YI <- ( outer(cI_2o_ap,p_XF*diffx) + outer(1-cI_2o_ap,p_XF_minus*diffx) ) 
  
  # extend this to 3d
  p_XF_cond_YIrep <- array(p_XF_cond_YI,c(resolution,resolution,resolution/2))
  # change order of dimensions for later integration (x: YI, y: reps, z: XF)
  p_XF_cond_YIrep <- aperm(p_XF_cond_YIrep,c(1,3,2))
  # Get the distribution of predicted XIs as a function of YI (x: YI, y: XI)
  p_XI_cond_YI_norm_t <- t(p_XI_cond_YI_norm) 
  # Get this XI distribution into 3d to allow for multiplication
  p_XI_cond_YI_rep <- array(p_XI_cond_YI_norm_t,c(resolution,resolution/2,resolution))
  # Multiply the two
  p_XIXF_cond_YI <- p_XI_cond_YI_rep*p_XF_cond_YIrep
  
  ## Computing the actual value
  # First, we get the array in the right orientation (x: XI, YI, XF)
  value_ZF_cubealt <- aperm(value_ZF_cube,c(2,1,3))  
  # We integrate over XI and XF (eq. 20)
  Q1 <- rowSums(p_XIXF_cond_YI*value_ZF_cubealt) + rS
  
  
  
  ## Statistics ---- ###########################################################
  
  
  ### Set-up ----
  
  ## Compute the Q-value difference
  Q_difference_n <- Q1 - Q0 # actual differences between the Q-values
  Q_difference <- sign(Q_difference_n) # binary difference
  diffcheck <- diff(Q_difference) 
  diffcheck <- c(diffcheck, diffcheck[length(diffcheck)-1]) # because diff cuts off one value
  
  
  # Get indices and values for when Q-values intersect
  index_threshold_YI <- which(diffcheck != 0)
  index_threshold_YI1 <- c(index_threshold_YI[1],
                           index_threshold_YI[1]+1)
  index_threshold_YI2 <- c(index_threshold_YI[2],
                           index_threshold_YI[2]+1)
  threshold_YI <- c(mean(x[index_threshold_YI1]),
                    mean(x[index_threshold_YI2]))
  
  
  ## Average seeking and confidence confidence thresholds ----
  
  if (rS == 0){ # if seeking is free
    
    av_search <- 1 # always seeking
    cI_cutoffs <- c(NA,NA)
    
  } else if (all(is.na(threshold_YI))) { 
    # if the Q-values never intersect
    
    if(min(Q_difference) > 0){
      av_search <- 1 # always seeking
    } else {
      av_search <- 0 # never seeking
    }
    # no cut-offs
    cI_cutoffs <- c(NA,NA)
    
    
    
  } else if (0 < max(Q_difference_n) && max(Q_difference_n) < Qdifftol || 0 > min(Q_difference_n) && min(Q_difference_n) > -Qdifftol ) {
    # if the maximal difference doesn't fall outside the the tolerance 
    # this takes effect when the Q-values values don't cross more than Qdifftol
    
    if (0 < max(Q_difference_n) && max(Q_difference_n) < Qdifftol) {
      av_search <- 0
    } else{
      av_search <- 1
    }
    
    cI_cutoffs <- c(NA,NA)
    
    
  } else {
    # in the normal case (when there is a reasonable intersection) >>>> 
    
    # First, the confidence thresholds
    if (is.na(threshold_YI[2])){ 
      # if we can only identify one cut-off (other one will lie outside of range )
      interval_YI_search <- c(-Inf,threshold_YI[1]) # let the interval run to -Inf
      cI_cutoffs <- c(NA,cI.2o(threshold_YI[1], sI_subjective, tI, rho,1)) # set up only one cI
      
    } else {
      # if we get two thresholds
      interval_YI_search <- threshold_YI
      cI_cutoffs <- cI.2o(interval_YI_search, sI_subjective, tI, rho, 1) 
      # note that this computes two thresholds due to the vectorization of cI.2o
      
    }
    
    # Compute the seeking averages by getting the cummulative probabliites
    p_search_plus <- pmvnorm(lower = c(0,interval_YI_search[1]), 
                             upper = c(Inf,interval_YI_search[2]), 
                             mean = c(d,d),
                             sigma = SI_objective)[[1]]
    p_search_minus <- pmvnorm(lower = c(-Inf,-interval_YI_search[2]), 
                              upper = c(0,-interval_YI_search[1]), 
                              mean = c(d,d),
                              sigma = SI_objective)[[1]]
    
    av_search <- p_search_plus + p_search_minus
    
  }
  
  
  ## Final accuracy ############################################################
  
  # NOTE:The subjective second-order is not optimized to support final accuracy
  # so this is commented out here
  
  #
  # if (aF == T) {
  #   
  #   if (av_search == 0) {
  #     # no final accuracy to compute when there is no seeking
  #     av_accuracy_s0 <- pnorm(0,-1,zetaI)
  #     av_accuracy_s1 <- NA
  #     
  #   } else {
  #     
  #     # Getting the ZFs: This allows us to check for a mistake
  #     # (per the sign of ZF)
  #     ZF_cube <- outer(ZI_matrix/zetaI^2, x/sF^2, FUN = "+")*zetaF^2
  #     # dim(ZF_cube)
  #     
  #     # The density of the individual values of ZF
  #     p_ZF_full <- outer(p_XIYI_matrix*diffxm,p_XF*diffx)
  #     # note that we now take the full XI space and not just half the space as before
  #     
  #     if (av_search == 1) { # if always seeking
  #       
  #       av_accuracy_s0 <- NA # no accuracy without seeking because there's always seeking
  #       
  #       # We still  need to compute the final accuracy
  #       ZF_cube <- outer(ZI_matrix/zetaI^2, x/sF^2, FUN = "+")*zetaF^2
  #       
  #       index_sought_correct <- ZF_cube > 0
  #       index_sought_incorrect <- ZF_cube < 0
  #       
  #       av_accuracy_s1 <- sum(index_sought_correct*p_ZF_full)
  #       incorrect_s1 <- sum(index_sought_incorrect*p_ZF_full)
  #       
  #       av_accuracy_s1 <- av_accuracy_s1 / (av_accuracy_s1 + incorrect_s1)
  #       # this could also be achieved by doing phi(zeta_F) but can
  #       # av_accuracy_s1 <- pnorm(0,-1,zetaF)
  #     }
  #     
  #     else {
  #       
  #       # Setting up indexes in XI - YI space where seeking happens
  #       index_sought_p <- outer(rep(1, length(diffcheck)/2), Q_difference)
  #       index_sought_m <- outer(rep(1, length(diffcheck)/2), rev(Q_difference))
  #       
  #       index_sought_full <- rbind(index_sought_m, index_sought_p)
  #       
  #       # set these indeces up in 3d
  #       index_sought_full_3d <- array(index_sought_full,c(resolution,
  #                                                         resolution,
  #                                                         resolution))
  #       
  #       # dim(p_ZF_full)
  #       # sum(p_ZF_full) # this should sum to ~ 1
  #       
  #       # Accuracy after seeking
  #       index_sought_correct <- ZF_cube > 0 & index_sought_full_3d == 1
  #       index_sought_incorrect <- ZF_cube < 0 & index_sought_full_3d == 1
  #       
  #       av_accuracy_s1 <- sum(index_sought_correct*p_ZF_full)
  #       incorrect_s1 <- sum(index_sought_incorrect*p_ZF_full)
  #       # normalisation
  #       av_accuracy_s1 <- av_accuracy_s1 / (av_accuracy_s1 + incorrect_s1)
  #       # av_accuracy_s1
  #       
  #       # Accuracy after not seeking (note how this is a happens in 2d (Z_I))
  #       index_not_sought_correct <- ZI_matrix > 0 & index_sought_full == -1
  #       index_not_sought_incorrect <- ZI_matrix < 0 & index_sought_full == -1
  #       # mass for inaccurate and accurate indeces
  #       av_accuracy_s0 <- sum(index_not_sought_correct*p_XIYI_matrix*diffxm)
  #       incorrect_s0 <- sum(index_not_sought_incorrect*p_XIYI_matrix*diffxm)
  #       # normalisation
  #       av_accuracy_s0 <- av_accuracy_s0 / (av_accuracy_s0 + incorrect_s0)
  #     }
  #   }
  # }
  
  
  
  ## Output ----- ##############################################################
  
  out <- list()
  if (output == "av_search") {
    out <- av_search
  } else if (output == "cI_cutoff") { 
    out <- cI_cutoffs
  } else if (output == "av_accuracy_s1") {
    out <- av_accuracy_s1
  } else if (output == "av_accuracy_s0") {
    out <- av_accuracy_s0
  } else if (output == "all"){
    out$av_search <- av_search
    out$cI_cutoffs <- cI_cutoffs
    if (aF == T) {
      out$av_accuracy_s1 <- av_accuracy_s1
      out$av_accuracy_s0 <- av_accuracy_s0
    }
  }
  
  return(out)
  
}

