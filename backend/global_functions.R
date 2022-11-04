# This defines or sources all the functions used throughout

# SOURCING -----

# set working directory (this usually doesn't need to be done, because it is
# called directly from the markdown)
# setwd("TheoryPaper/paper_code/") 


# Modelling
source("backend/models/joint/model_functions_joint.R")
source("backend/models/postdecisional/model_postdecisional.R")
source("backend/models/secondorder/model_2o.R")
source("backend/models/secondorder/model_2o_functions.R")

# Plotting
source("backend/plotting/plot_av_seek_heatmap.R")
source("backend/plotting/plot_metacognitive_efficiency.R")
source("backend/plotting/plot_confidence_distributions.R")



# Global functions ----
# Functions that are used throughout the documents

# Normalizes a vector or matrix to sum to 1
normalize <- function(object){
  return(object/sum(object))
}


# Progress bar functions
qpb <- function(df){ # set-up progressbar based on dataframe (length = nrow)
  txtProgressBar(min = 0, max = nrow(df), initial = 0, char = "//", style = 3) 
}

upb <- function(pb,i){ # update that progress bar
  setTxtProgressBar(pb,i)}


# Function that does expand.grid but adds Colnames functionality
expand.grid.df <- function(...,col_names = F){
  
  df <- expand.grid(...)
  
  if (col_names[1] != F) {
    colnames(df) <- col_names
  }
  return(df)
}


`%!in%` <- Negate(`%in%`)
