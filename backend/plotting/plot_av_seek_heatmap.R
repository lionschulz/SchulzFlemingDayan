plot_av_seek_heatmap <- function(dataframe, tI = T, legend = F, size = "big"){
  
  plot <- ggplot(dataframe, aes(x = phi_sI, y = phi_sF)) + 
    geom_tile(aes(fill=av_search)) +
    scale_x_continuous(name=bquote("Initial Accuracy. " * phi *"("*sigma["I"]*")"),
                       expand = c(0,0)) +
    scale_y_continuous(name=bquote(X[F] * " Accuracy " * phi *"("*sigma["F"]*")"),
                       expand = c(0,0)) +
    scale_fill_viridis_c(option = "magma", limits = c(0,1),
                         breaks = c(0,.5,1),
                         name = "Av.\nSearch") 
  
  if (legend == F) plot <- plot + theme(legend.position = "none")
  
  if (tI == T) tI_temp <- dataframe$tI[1]
  if (tI == T) plot <- plot + ggtitle(bquote(tau["I"] == .(tI_temp)))
  
  if (size == "small") {
    
    breaks_small <- c(.6,.9)
    labels_small <- c(".6",".9")
    
    plot <- plot + 
      scale_y_continuous(name=bquote(X[F] * " Accuracy " * phi *"("*sigma["F"]*")"),
                         expand = c(0,0),
                         breaks = breaks_small,
                         labels = labels_small) +
      scale_x_continuous(name=bquote("Initial Accuracy " * phi *"("*sigma["I"]*")"),
                         expand = c(0,0),
                         breaks = breaks_small,
                         labels = labels_small) 
    
  }
  return(plot)
}
