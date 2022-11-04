# Functions that are used when approaching the models through a simulation lens


# Returns a df of simulated random XI and YI given the model parameters
simulate_XIYI <- function(
  sI,
  tI,
  rhoI = 0, # 0 for the postdecisional model
  n_total_simulations = 1000
){

  # print(n_total_simulations)
  
  # 1st - we set up the parameters of the stimuli
  SI <- computeCov(sI, tI, rhoI)
  d_vector <- c(1, 1)  
  
  # 2nd - we simulate a number of random variables from this parameter set
  simulated_XIYI_positive <- rmvnorm(n = n_total_simulations, 
                                     mean = d_vector, 
                                     sigma = SI)
  simulated_XIYI_negative <- rmvnorm(n = n_total_simulations, 
                                     mean = -1 * d_vector, 
                                     sigma = SI)
  all_XIYI <- rbind(simulated_XIYI_positive, simulated_XIYI_negative)
  
  # 3rd - we get these random variables into a dataframe ...
  df_simulations <- data.frame(d = rep(c(1, -1), each = n_total_simulations),
                               XI = all_XIYI[, 1], 
                               YI = all_XIYI[, 2])
  
  return(df_simulations)
}


# Computes the action and confidence based on a dataframe of XI's and YI's 
# (produced by 'simulate_XIYI')
get_aI_and_cI_from_df_simulations <- function(
  df_simulations,
  sI, # note that this sigma doesn't need to be the one in simulate_XIYI <- over/underconfidence
  tI,
  rhoI = 0, # 0 for the postdecisional model
  model = "2o", # "Pd" or "2o"
  sI_subjective = FALSE
) {
  
  # ... and then compute the decision and confidence based on this
  df_simulations$aI <- sign(df_simulations$XI)
  df_simulations$correct <- df_simulations$aI == df_simulations$d
  
  # confidence depending on the model
  if (model == "2o") {
    df_simulations$cI <- cI.2o(df_simulations$YI,
                               sI,tI,rhoI,
                               df_simulations$aI)
  } else if (model == "Pd") {
    zetaI2 <- joinsd(sI, tI)^2 # combined variance
    df_simulations$zi <- (df_simulations$XI/sI^2 + df_simulations$YI/tI^2) * zetaI2 # combined cue 
    df_simulations$cI <- sigmoid(2 * df_simulations$aI * df_simulations$zi / zetaI2) # get the confidence conditioned on ai (i.e. XI)
  }
  
  return(df_simulations)
}



# Prepare df_simulations for the meta-d' calculations
prep_df_sim_for_metad <- function(
  df_simulations, 
  pad_factor, 
  padded
){
  
  # We define the borde of the confidence space
  confidence_bin_border <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  n_ratings <- length(confidence_bin_border)
  
  # We apply these borde to the data, fit in a way that is
  # not conditioned on choice
  df_simulations$cI_bin <- .bincode(df_simulations$cI,
                                    breaks = confidence_bin_border)
  # ... and then in a way that is conditioned on choice
  df_simulations$cI_direction <- df_simulations$cI_bin * df_simulations$aI
  
  # 7th we get the data using the format used by the meta-d' function
  df_from_counts <- suppressMessages(
    sdt_counts(df_simulations, stimulus = d, response = cI_direction)
  )
  # Confidence ratings and not all - so we need to fill the empty cells
  df_counts_base <- data.frame(cI_direction = c(-5:-1, 1:5))
  df_counts_base$"-1" <- 0
  df_counts_base$"1" <- 0
  # The problem with sdt_counts is that it only counts the empirically existing 
  # confidence ratings and not all - so we need to fill the empty cells
  df_counts_base <- data.frame(cI_direction = c(-5:-1, 1:5))
  df_counts_base$"-1" <- 0
  df_counts_base$"1" <- 0
  df_counts_base
  
  df_counts_base[df_counts_base$cI_direction %in% df_from_counts$cI_direction, 
  ]$"-1" <- df_from_counts$"-1"
  df_counts_base[df_counts_base$cI_direction %in% df_from_counts$cI_direction,
  ]$"1" <- df_from_counts$"1"
  # We do some padding to facilitate fitting
  df_counts_base$"-1" <- df_counts_base$"-1" + pad_factor
  df_counts_base$"1" <- df_counts_base$"1" + pad_factor
  
  
  
  # We get the vectors in the right format
  if (padded) {
    nR_S1 <- df_counts_base$"-1"
    nR_S2 <- df_counts_base$"1" 
  } else {
    nR_S1 <- df_from_counts$"-1"
    nR_S2 <- df_from_counts$"1" 
  }

  # print(nR_S1)
  # print(nR_S2)
 
  return(
    list(nR_S1 = nR_S1,
         nR_S2 = nR_S2)
  )
}

# Run the meta-d' prime analysis on a df_simulations dataframe produced by 
# simulate_XIYI() and processed by get_aI_and_cI_from_df_simulations()
get_metad_from_df_simulations <- function(
  df_simulations, 
  pad_factor, 
  padded,
  fit_type = "MLE",
  hmetad_seed = FALSE
) {
  
  vectors_for_meta_d <- prep_df_sim_for_metad(
    df_simulations, 
    pad_factor, 
    padded
  )
  
  # print(vectors_for_meta_d)
  
  if(fit_type == "MLE"){
    # fit meta-d'
    fit_MLE <- fit_meta_d_MLE(as.numeric(vectors_for_meta_d$nR_S1), 
                              as.numeric(vectors_for_meta_d$nR_S2))
    
    fit <- list()
    fit$meta_dprime <- fit_MLE$meta_da[1]
    fit$dprime <- fit_MLE$da[1]
    fit$mratio <- fit_MLE$M_ratio[1]
    
  } else if(fit_type == "hmetad"){
    
    fit <- fit_h_meta_d(as.numeric(vectors_for_meta_d$nR_S1), 
                        as.numeric(vectors_for_meta_d$nR_S2),
                        seed = hmetad_seed)
    
  }

  return(
    list(fit = fit, data = vectors_for_meta_d)
  )
}



# Function that gets meta-d' and related parameters for a simulated confidence
# model 
# Parameters:  
# sI, tI, rhoI are parameters of the model
# n_total_simulations defines the number of simulations to run
# pad_factor pads the confidence vectors to facilitate fitting
# paddded controls whether to use the padded vector or the 'empirical' vector
# just using the confidence judgements used in the data
get_meta_parameters_from_model_parameters <- function(
  sI,
  tI,
  rhoI = 0, # 0 for the postdecisional model
  n_total_simulations = 1000,
  pad_factor = 1,
  padded = TRUE,
  model, # "Pd" or "2o",
  fit_type = "MLE"
) {
  
  # Make sure that rhoI == 0 for postdecisional model
  rhoI <- ifelse(model == "Pd", 0, rhoI)
  
  # Simulate the random variables
  df_simulations <- simulate_XIYI(sI, tI, rhoI, n_total_simulations)
  
  # Get the confidence and action
  df_simulations <- get_aI_and_cI_from_df_simulations(df_simulations, 
                                                      sI, tI, rhoI, 
                                                      model)
  
  metad_results <- get_metad_from_df_simulations(
    df_simulations, 
    pad_factor, 
    padded,
    fit_type
    )
  return(metad_results)
}



# Gets metacognitive and search statistics from a simulated model
# Works for both postdecisional and second-order model
# Input: model parameters, model type, potential subjective sigma_I (for 
# second order model only)
# n_total_simulations defines the number of trials to simulate
# Ouput: list with results from the numerical modelling and the 'on-policy' simulations
get_search_and_meta_stats <- function(
  sI, tI, rhoI, 
  sF, rS,
  model = "2o",
  sI_subjective = FALSE,
  n_total_simulations = 1000,
  pad_factor = 1, # for fitting meta-d',
  compute_meta_d = TRUE,
  meta_d_fit_type = "MLE",
  hmetad_seed = FALSE,
  compute_search_stats = TRUE
){
  
  # First-order statistic
  av_correct_expected <- pnorm(0, -1, sI)

  # Simulate initial action and confidence
  rhoI <- ifelse(model == "2o", rhoI, 0)
  df_simulations <- simulate_XIYI(sI, tI, rhoI, n_total_simulations)
  
  sI_for_confidence <- ifelse(sI_subjective, sI_subjective, sI)
  df_simulations <- get_aI_and_cI_from_df_simulations(df_simulations, 
                                                      sI_for_confidence, tI, rhoI, 
                                                      model, 
                                                      sI_subjective)
  # Compute summary statistics for this
  av_correct_empirical <- mean(df_simulations$correct)
  av_confidence_empirical <- mean(df_simulations$cI)
  
  if(compute_meta_d){
    metad_results <- get_metad_from_df_simulations(df_simulations, 
                                                   pad_factor = pad_factor,
                                                   padded = TRUE,
                                                   fit_type = meta_d_fit_type,
                                                   hmetad_seed = hmetad_seed)
  }
  
  
  # Search statistics
  if(compute_search_stats){
    
    if (model == "2o") {
      if (sI_subjective) {
        search_stats <- model.2o_subjective(
          sI_objective = sI,
          sI_subjective = sI_subjective,
          tI = tI,
          rho = rhoI,
          sF = sF,
          rS = rS,
          subjective_confidence = TRUE,
          resolution = 180,
          output = "all"
        )      
      } else {
        search_stats <- model.2o(sI, tI, rhoI, sF, rS, 
                                 output = "all",resolution = 150)
      }
      
      search_stats$av_search_when_correct <- search_stats$av_search_when_correct / av_correct_expected
      search_stats$av_search_when_incorrect <- search_stats$av_search_when_incorrect / (1 - av_correct_expected)
      
    } else if (model == "Pd") {
      zetaI <- joinsd(sI, tI)
      search_stats <- model.Pd(zetaI, sF, rS, aF = F, resolution = 1000)
    }
    
    # Compute whether agent seeks on a trial by trial basis
    # use only upper cut-off for better symmetry in second-order model
    search_stats$cI_cutoff <- ifelse(model == "2o", search_stats$cI_cutoffs[2], search_stats$cI_cutoff)
    df_simulations$search <- df_simulations$cI > (1 - search_stats$cI_cutoff) & 
      df_simulations$cI < search_stats$cI_cutoff
    
    
    # Investigate the conditional search
    empirical_search_by_correct <- df_simulations %>%
      group_by(correct) %>%
      summarise(mean = mean(search))
    empirical_av_search_correct <- empirical_search_by_correct$mean[2]
    empirical_av_search_incorrect <- empirical_search_by_correct$mean[1]
    
    empirical_search <- list(
      average = mean(df_simulations$search),
      when_correct = empirical_av_search_correct,
      when_incorrect = empirical_av_search_incorrect
    )
    
    # Compute the search slope (Seek ~ Accuracy on trial by trial basis)
    # glm_summary <- summary(
    #   bayesglm(search ~ as.integer(correct), 
    #            data = df_simulations,
    #            family = "binomial")
    # )
  }

  # Get output ready
  output_list <- list(
    av_correct = list("expected" = av_correct_expected,
                      "empirical" = av_correct_empirical),
    av_confidence = av_confidence_empirical,
    # glm_summary = glm_summary,
    df_simulations = df_simulations
  )
  
  if(compute_search_stats){
    output_list$search_stats <- search_stats
    output_list$empirical_search <- empirical_search
  }
  
  # Compute meta-d'
  if(compute_meta_d){
    output_list$metad_results <- metad_results
  }

  return(output_list)
}

# Computes the search sensitivity and other search/metacog measures based on a
# vector of sigma_Is and tau_Is (and potentially subjective sigma_Is). Works
# for both the second-order and postdecisional model
get_simulated_df_with_search_and_meta_metrics <- 
  function(sI_vector, 
           tI_vector,          
           sI_subj_vector = FALSE,
           rho_fix, 
           sF_fix, 
           rS_fix, 
           model = "2o",
           n_total_simulations = 10000,
           compute_meta_d = TRUE,
           pad_factor = 3,
           meta_d_fit_type = "MLE",
           hmetad_seed = FALSE){

  # set up the dataframe ----
  df <- data.frame(
    id = 1:length(sI_vector),
    sI = sI_vector,
    sI_subj = sI_subj_vector,
    tI = tI_vector,
    sI_by_tI = sI_vector / tI_vector,
    av_correct = NA, 
    av_confidence = NA,
    av_search = NA,
    av_search_when_correct = NA,
    av_search_when_incorrect = NA,
    da = NA,
    meta_da = NA,
    M_ratio = NA,
    search_slope = NA,
    av_search_empirical = NA,
    av_search_correct_empirical = NA,
    av_search_incorrect_empirical = NA
  )
  
  # check if we use a subjective sigma_I (for over-/underconfidence) ----
  if(sI_subj_vector[1]){
    df <- df %>% 
      mutate(
        sigma_ratio = sI / sI_subj,
        sigma_difference = sI - sI_subj
      )
  }
  
  # get the actual metrics ----
  for (i in 1:nrow(df)) {
    
    cat(paste(i, "/", nrow(df), " | "))
    
    metric_list <- get_search_and_meta_stats(
      sI = df$sI[i],
      tI = df$tI[i],
      rhoI = rho_fix,
      sF= sF_fix,
      rS = rS_fix,
      model = model,
      sI_subjective = df$sI_subj[i],
      n_total_simulations = n_total_simulations, 
      pad_factor = 2,  # note that this is quite high but helps to work with the 
                       # extreme parameter combinations we're providing
      compute_meta_d = compute_meta_d,
      meta_d_fit_type = meta_d_fit_type,
      hmetad_seed = hmetad_seed
    )
    
    # print(metric_list)
    
    # add the metrics to the dataframe ----
    df$expected_accuracy <- pnorm(0, 1, sd = sI_vector)
    df$av_correct[i] <- metric_list$av_correct$empirical
    df$av_confidence[i] <- metric_list$av_confidence 
    df$av_search[i] <- metric_list$search_stats$av_search
    if(model == "2o" & all(sI_subj_vector == FALSE)){
      # only use the numerical method to get conditional searches when
      # the model's confidence is calibrated
      df$av_search_when_correct[i] <-
        metric_list$search_stats$av_search_when_correct
      df$av_search_when_incorrect[i] <-
        metric_list$search_stats$av_search_when_incorrect
    }

    if(compute_meta_d){

      df$da[i] <- metric_list$metad_results$fit$dprime
      df$meta_da[i] <- metric_list$metad_results$fit$meta_dprime
      df$M_ratio[i] <- metric_list$metad_results$fit$mratio
      # df$search_slope[i] <- metric_list$glm_summary$coefficients[2,1]
      
    }

    df$av_search_empirical[i] <- metric_list$empirical_search$average
    df$av_search_correct_empirical[i] <- metric_list$empirical_search$when_correct 
    df$av_search_incorrect_empirical[i] <- metric_list$empirical_search$when_incorrect 
    

  }
  
  df <- df %>% 
    mutate(search_ratio = 
             (av_search_when_incorrect - av_search_when_correct) / av_search,
           search_ratio_empirical = 
             (av_search_incorrect_empirical - av_search_correct_empirical) /
             av_search_empirical,
           search_ratio_empirical2 = 
             (av_search_incorrect_empirical - av_search_correct_empirical) /
             av_search,
           search_slope = -search_slope)
  
  return(df)
}



# Plots a 'histogram' of the confidence ratings of the data provided to the
# fit_meta_d_MLE() functions
plot_confidence_for_metad <- function(nR_S1, nR_S2, facetted = FALSE){
  
  n_ratings <- length(nR_S1)/2
  
  df_plot_meta <- data.frame(
    stim = rep(c("1", "2"), each = length(nR_S1)),
    response = rep(rep(1:2,each = n_ratings),2),
    rating = rep(c(n_ratings:1,1:n_ratings),2),
    correct = NA,
    count = c(nR_S1, nR_S2)
  )
  
  df_plot_meta$correct <- df_plot_meta$stim == df_plot_meta$response
  
  pl <- df_plot_meta %>% 
    ggplot(aes(x = rating, y = count, fill = correct)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c(colour$distribution_incorrect,
                                   colour$distribution_correct))
    
    if (facetted) {
      return(pl  + facet_grid(stim~.))
    } else {
      return(pl)
    }
}
