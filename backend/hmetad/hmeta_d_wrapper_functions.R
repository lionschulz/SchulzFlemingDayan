# Wrapper function for the individual metad fit function from h-meta-d'
fit_h_meta_d <- function(nR_S1, nR_S2, seed = FALSE){
  
  # Apply the function
  fit_h <- metad_indiv(nR_S1, nR_S2, seed)
  
  print("fit done")
  
  # Mean values 
  Value <- summary(fit_h$samples)
  stat <- data.frame(mean = Value[["statistics"]][, "Mean"])
  stat %<>%
    rownames_to_column(var = "name")
  
  
  meta_d <- stat[which(stat$name == "meta_d"),]$mean
  dprime <- fit_h$data$d1
  mratio <- meta_d / dprime
  
  results = list(
    meta_dprime = meta_d,
    dprime = dprime,
    mratio = mratio
  )
  
  return(results)
}