# This file contains the functions that are used in both models

# Function that takes two standard deviation and returns their joint sd
joinsd <- function(sd1, sd2){
  sdout <- 1 / ( 1 / sd1^2 + (1 / sd2^2))
  sdout <- sqrt(sdout)
}

# Function used to compute posteriors/confidences in the models
sigmoid <- function(x){
  1/(1+exp(-x))
}