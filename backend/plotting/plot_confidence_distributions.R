# this function produces the confidence distribution figures
# both for the postdecisional and the second-order moodels
# 
# sI = 2
# tI = 2
# rho = .3
# sF = 1
# rS = -.1
# showrho = TRUE
# xlabel = TRUE
# size_settings = TRUE
# model = "2o"
# show_title = TRUE

plot_cI_distribution <- function(model, sI,tI,rho = 0,sF = FALSE,rS = FALSE,showrho = FALSE, xlabel = FALSE, size_settings = TRUE, show_title = TRUE, sI_subj = FALSE){

  # to create the figures we simulate nsim random variables from the distribution
  # and then plot the respective confidences
  
  nsim <- 300000 # number of simulations
  d <- c(1,1) # mean of source (only use first for the postdecisional model)
  
  
  # Simulating the confidence distributions ----
  
  if (model == "Pd") {
    # simulate the random variable values
    XI_vector <- rnorm(nsim, d, sI)
    
    if (tI < Inf) { # in case we want the first-order model
      YI_vector <- rnorm(nsim, d, tI)
    } else { # because rnorm doesn't support sd = Inf
      YI_vector <- rep(0,nsim)
    }
    
    # get the confidence
    zetaI2 <- joinsd(sI, tI)^2 # combined variance
    ZI <- (XI_vector/sI^2 + YI_vector/tI^2)*zetaI2 # combined cue 
    cI_vector <- sigmoid(2*sign(XI_vector)*ZI/zetaI2) # get the confidence conditioned on ai (i.e. xi)
  
  } else if (model == "2o"){
  
    # covariance matrix
    Si <- computeCov(sI,tI,rho)
    
    # simulate YIs and XIs
    XIYI <- rmvnorm(n = nsim, mean = d, sigma = Si) # samples from bivariate normal
    XI_vector <- XIYI[,1] 
    YI_vector <- XIYI[,2]
    
    # compute the confidence
    aI <- sign(XI_vector)
    
    if(sI_subj == FALSE){
      cI_vector <- cI.2o(XIYI[,2],sI,tI,rho,aI)
    } else {
      cI_vector <- cI.2o(XIYI[,2],sI_subj,tI,rho,aI)
    }
 
    
  }
  
  # Getting the seeking zones (if wanted, else: sF = rS = FALSE) ---
  if (sF != FALSE) {
    if (model == "Pd") {
      # setting up the search zone and the connfidence cut-offs
      zetaI <- sqrt(zetaI2)
      cI_cutoff <- model.Pd(zetaI,
                            sF,
                            rS,
                            aF = FALSE)$cI_cutoff
      
      
      if (tI < Inf) { # in normal postdecisional models, we mirror the upper threshold
        cI_cutoff_low <- 1 - cI_cutoff
      } else { # the lower threshold is at 50 % confidence for the first-order model
        cI_cutoff_low <- 0.5
      }
    } else if (model == "2o"){
      
      cI_cutoff <- model.2o(sI,tI,rho,sF,rS,output = "cI_cutoff")
      
      # move the line out of the way when there is no probing cut-off
      if (all(is.na(cI_cutoff))) cI_cutoff <- 100 
      # make opposite line when out of bounds
      if (is.na(cI_cutoff[1])) cI_cutoff[1] <- 1 - cI_cutoff[2]
      
      # adapt variable names to conform with postdecisional plot
      cI_cutoff_low <- cI_cutoff[1]
      cI_cutoff <- cI_cutoff[2]
    }
  }
  
  # Getting this into a dataframe ----
  
  df_cIs <- data.frame(XI = XI_vector, YI = YI_vector, cI = cI_vector)
  # we mark the action as correct and incorrect through their XIs
  df_cIs$aI <- "correct"
  df_cIs[df_cIs$XI < 0,]$aI <- "incorrect" 
  
  
  # Do the actual plotting ----
  
  # this allows us to adjust the height of the annotation box denoting
  # the seeking zone properly
  binwidth <- .025
  histbreaks <- seq(0,1,by = binwidth)
  histpre <- hist(df_cIs[sign(df_cIs$XI) == 1,]$cI, plot = FALSE, breaks = histbreaks)
  plotheight <- max(histpre$counts)
  
  # base plot only with the distribution
  plot <- ggplot(df_cIs, aes(x = cI, fill = aI)) + 
    geom_histogram(breaks = histbreaks,position = "identity", alpha = .7) +
    scale_fill_manual(values = c(colour$distribution_correct,colour$distribution_incorrect)) + 
    xlim(c(-.025,1.025)) +
    theme(axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y =element_blank()) +
    guides(fill = "none") +
    labs(y = "", x = "Confidence") +
    scale_y_continuous(expand=c(0,0)) 
  
  
  # Plotting settings ----
  
  if (sF != FALSE) {
    plot <- plot + 
      annotate(geom="rect", xmin =cI_cutoff_low, xmax = cI_cutoff, ymin = 0, 
               ymax = plotheight, fill = colour$distribution_area, alpha = .1) +
      geom_vline(xintercept = c(cI_cutoff_low, cI_cutoff), 
               colour = colour$distribution_area, 
               size = 2,
               alpha = colour$distribution_line_alpha) 
  }
  
  # if I want to print the parameter settings:
  if (show_title == TRUE) {
    if (showrho == TRUE) {
      title <- bquote(sigma["I"] == .(sI) ~ ", " ~
                        tau["I"] == .(tI) ~ ", " ~
                        rho["I"] == .(rho))
    } else {
      title <- bquote(sigma["I"] == .(sI) ~ ", " ~
                        tau["I"] == .(tI))
    }
    
    plot <- plot + ggtitle(title) 
  }
  
  if (size_settings == TRUE) {
    plot <- plot & scale_x_continuous(limits = c(0,1),
                                      breaks = c(0,.25,.5,.75,1),
                                      labels = c("0","","Confidence","","1")) &
      theme(plot.title = element_text(size = .8*base_size))
  }
  
  if (xlabel == FALSE) {
    plot <- plot + theme(axis.title.x=element_blank())
  }
  
  
  
  return(plot)

}










# 
# # second-order ----
# 
# 
# 
# # simulate YIs and XIs
# XIYI <- rmvnorm(n = nsim, mean = d, sigma = Si)
# # hist(XIYI[,1])
# # hist(XIYI[,2])
# ai <- sign(XIYI[,1])
# 
# # compute the confidence
# cI <- cImeta(XIYI[,2],sI,tI,rho,ai)
# # hist(cI)
# # mean(cI)
# 
# # organize everything
# df_dis <- data.frame(XI = XIYI[,1], YI = XIYI[,2], cI)
# 
# df_dis$aI <- "correct"
# df_dis[df_dis$XI < 0,]$aI <- "incorrect"
# 
# # settings for the histogram
# binwidth <- .025
# histbreaks <- seq(0,1,by = binwidth)
# histpre <- hist(df_dis[sign(df_dis$XI) == 1,]$cI, plot = FALSE, breaks = histbreaks)
# plotheight <- max(histpre$counts)*1.1
# 
# # make plot
# plot <- ggplot(df_dis, aes(x = cI, fill = aI)) + 
#   geom_histogram(breaks = histbreaks,position = "identity", alpha = .7) +
#   scale_fill_manual(values = c(colour$distribution_correct,colour$distribution_incorrect)) + 
#   xlim(c(-.025,1.025)) +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.line.y =element_blank()) +
#   guides(fill = "none") +
#   labs(y = "", x = "Confidence") +
#   scale_y_continuous(expand=c(0,0)) +
#   
#   
#   if (size_settings == T) {
#     plot <- plot & scale_x_continuous(limits = c(0,1),
#                                       breaks = c(0,.25,.5,.75,1),
#                                       labels = c("0","","Confidence","","1")) &
#       theme(plot.title = element_text(size = .8*base_size))
#   }
# 
# 
# 
# # if I want to plot the confidence cut-off
# if (sF != F) {
#   cI_cutoff <- seekPercent_2o(sI,tI,rho,sF,rS,output = "cIcut")
#   
#   # move the line out of the way when there is no probing cut-off
#   if (all(is.na(cI_cutoff))) cI_cutoff <- 100 
#   # make opposite line when out of bounds
#   if (is.na(cI_cutoff[1])) cI_cutoff[1] <- 1 - cI_cutoff[2]
#   
#   
#   
#   # add this information to the plot
#   plot <- plot + 
#     geom_vline(xintercept = cI_cutoff,
#                colour = colour$distribution_area, size = 2,
#                alpha = colour$distribution_line_alpha) +
#     annotate(geom="rect", xmin = cI_cutoff[1], xmax = cI_cutoff[2], ymin = 0, 
#              ymax = plotheight, fill = colour$distribution_area, alpha = .1)
# }
# 
# # if I want to print the parameter settings:
# if (show_title == T) {
#   if (showrho == T) {
#     title <- bquote(sigma["I"] == .(sI) ~ ", " ~
#                       tau["I"] == .(tI) ~ ", " ~
#                       rho["I"] == .(rho))
#   } else {
#     title <- bquote(sigma["I"] == .(sI) ~ ", " ~
#                       tau["I"] == .(tI))
#   }
#   
#   plot <- plot + ggtitle(title) 
#   
# }




# plot <- ggplot(df_cIs, aes(x = cI, fill = aI)) + 
#   geom_histogram(breaks = histbreaks,position = "identity", alpha = .7) +
#   scale_fill_manual(values = c(colour$distribution_correct, colour$distribution_incorrect)) +
#   geom_vline(xintercept = c(cI_cutoff_low, cI_cutoff), 
#              colour = colour$distribution_area, 
#              size = 2,
#              alpha = colour$distribution_line_alpha) +
#   xlim(c(-.025,1.025)) +
#   guides(fill = FALSE) +
#   labs(y = "") +  scale_y_continuous(expand=c(0,0)) +
#   theme(axis.text.y=element_blank(),
#         axis.title = element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.line.y =element_blank()) +
#   annotate(geom="rect", xmin =cI_cutoff_low, xmax = cI_cutoff, ymin = 0, 
#            ymax = plotheight, fill = colour$distribution_area, alpha = .1) +
#   ggtitle(bquote(sigma["I"] == .(sI) ~ ", " ~ tau["I"] == .(tI)))