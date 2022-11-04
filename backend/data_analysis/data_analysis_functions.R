get_averages_and_meta_stats <- function(df, df_add, fit_type = "MLE", hmetad_seed = FALSE){
  # Computes the averages from the main df as well as the relevant meta
  # statistics. Merges this into one summary df with the additional information
  # provided by df_add
  
  # Compute averages and merge with df_add ----
  df_averages <- df %>% 
    mutate(Confidence1st = as.numeric(Confidence1st)) %>% 
    group_by(ID) %>% 
    summarise(mean_accuracy = mean(Correct1st),
              mean_search = mean(SeeAgain),
              mean_confidence = mean(Confidence1st),
              mean_final_accuracy = mean(Correct2nd))
  
  df_averages <- merge(df_averages, df_add_complete)
  
  
  # Get the data ready for meta-d' ----
  
  ids <- unique(df_complete$ID)
  
  # labelling the side picked 
  df %<>% 
    mutate(initial_action = ((Correct1st*2) - 1) * Target) 
  
  df$initial_confidence <- df$initial_action * as.numeric(df$Confidence1st)
  
  df_confidence_count_data <- df %>% 
    group_by(ID, Target,initial_confidence) %>% 
    count()
  df_confidence_count_data
  
  df_confidence_count_setup <- data.frame(
    ID = rep(ids, each = 12), 
    Target = rep(c(-1, 1), each = 6),
    initial_confidence = c(-3:-1,1:3))
  df_confidence_count_setup
  
  df_confidence_count <- merge(df_confidence_count_setup, 
                               df_confidence_count_data, all = TRUE)
  df_confidence_count[is.na(df_confidence_count$n),]$n <- 0
  df_confidence_count
  
  df_for_meta <- df_confidence_count %>% 
    arrange(ID) %>% 
    select(ID, Target, initial_confidence, n) %>% 
    spread(initial_confidence, n)
  df_for_meta
  
  
  # Computing meta-d' ----
  df_metas <- data.frame(ID = ids, dprime = NA, meta_dprime = NA, mratio = NA)
  
  for (i in 1:length(ids)) {
    
    # grab id and data
    id_temp <- ids[i]
    df_meta_for_one_participant <- df_for_meta %>% 
      filter(ID == id_temp)
    
    if(fit_type == "MLE"){
      # fit meta-d'
      fit_MLE <- fit_meta_d_MLE(as.numeric(df_meta_for_one_participant[1,3:8]), 
                                as.numeric(df_meta_for_one_participant[2,3:8]))
      # fit_MLE
      
      df_metas[i,]$meta_dprime <- fit_MLE$meta_da[1]
      df_metas[i,]$dprime <- fit_MLE$da[1]
      df_metas[i,]$mratio <- fit_MLE$M_ratio[1]
      
    } else if(fit_type == "hmetad"){
      
      fit_hierarchical <- fit_h_meta_d(as.numeric(df_meta_for_one_participant[1,3:8]), 
                                       as.numeric(df_meta_for_one_participant[2,3:8]),
                                       hmetad_seed)
      
      df_metas[i,]$meta_dprime <- fit_hierarchical$meta_d
      df_metas[i,]$dprime <- fit_hierarchical$dprime
      df_metas[i,]$mratio <- fit_hierarchical$mratio
      
    }
    
    cat(paste(i, "/", length(ids), "|"))
  }

  # Output ----
  df_summary <- merge(df_averages, df_metas)
  return(df_summary)
}