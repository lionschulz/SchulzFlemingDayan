# This is the function that plots the average search by metacognitive efficiency (by accuracy plot)

plot_meta_efficiency <- function(df, model){
  
  if (model == "Pd") {
    xaxis_label <- bquote(sigma[I]*"/"*zeta[I])
  } else if (model == "2o"){
    xaxis_label <- bquote(sigma[I]*"/"*tau[I])
  }
  
  plot <- ggplot(df, aes(x = efficiency_ratio, av_search, colour = as.factor(sI))) +
    geom_vline(xintercept = 1, size =lsize, colour = "grey") +
    geom_line(size = lsize) +
    geom_point(size = psize) +
    ylim(c(0,1)) +
    scale_color_manual(name = bquote(sigma[I]), 
                       values = palette_for_sigma) + 
    labs(y = "Av. Search", x = xaxis_label) +
    scale_x_continuous(limits = c(.6,1.8), breaks = c(.6,1,1.4,1.8)) +
    scale_y_continuous(limits = c(0,1),breaks = c(0,.5,1))
  plot
}
