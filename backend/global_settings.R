# File that defines the settings (mostly for plotting) used across the markdowns

# Rendering ----
# This enables anti-aliasing on Windows

# Enable anti-aliasing on Windows
if(Sys.info()['sysname'] == "Windows"){
  
  trace(grDevices::png, quote({
    if (missing(type) && missing(antialias)) {
      type <- "cairo-png"
      antialias <- "subpixel"
    }
  }), print = FALSE)
  
  
  # Enable anti-aliasing on Windows
  trace(grDevices:::png, quote({
    if (missing(type) && missing(antialias)) {
      type <- "cairo-png"
      antialias <- "subpixel"
    }
  }), print = FALSE)
  
  
}

# A function that joins together lines and points, so that 
# we can avoid weird overlap issues
# taken from: 
# https://stackoverflow.com/questions/42983499/plot-points-in-front-of-lines-for-each-group-ggplot2-equivalent-of-type-o
linepoint = function(data, group.var, lsize=1.2, psize=4) {
  lapply(split(data, data[,group.var]), function(dg) {
    list(geom_line(data=dg, size=lsize),
         geom_point(data=dg, size=psize)) 
  })
}


# Basics ----
base_size <- 30
theme_set(theme_classic(base_size = base_size))


# Sizes ----
## Results
lsize <- 1.5
psize <- 3

## Introduction
lsize_intro <- 2
size_annotate <- 10

# Colours ----

# Single colours ----
colour <- list()
colour$av_confidence <- "darkgrey" # defines colour of confidence line
colour$Q0 <- "grey"
colour$posterior_ZI <- "grey45"

colour$distribution_correct <-"chartreuse2"
colour$distribution_incorrect <- "red2"
colour$distribution_area <- "black"
colour$distribution_line_alpha <- .5

# Palettes -----

palette_5taus <- rev(brewer.pal(n = 8, name = "Blues"))
palette_4costs <- brewer.pal(n = 5, name = "YlGn")[2:5]
palette_for_sigma <- sequential_hcl(5, "SunsetDark")
palette_s1_corr <- rev(brewer.pal(6,"YlOrRd")[2:6]) 
palette_2ais <- brewer.pal(6,"Greys")[c(5,3)]



# Shapes
shape_s1_corr <- 15


# Axes ---- 
breaks_avseek_y <- c(0,.5,1)
axis_name_initial_acc <- bquote("Initial Acc. " * phi *"("*sigma["I"]*")")

