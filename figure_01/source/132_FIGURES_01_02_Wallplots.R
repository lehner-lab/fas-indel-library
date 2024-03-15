library(tidyverse)
library(scales)
library(pheatmap)
load("data/indels_library.RData")



idx_singles <- which(nchar(indels_library$sequence)==63 & indels_library$id != "WT")
density_object <- density(x = indels_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]
# plot(density_object)
# abline(v=library_bias)
# plot(density(indels_library$es - library_bias))
# abline(v=0)
# abline(v=1)
indels_library$es <- indels_library$es - library_bias

indels_library$psi <- (49.1/1) * exp(indels_library$es)
indels_library$psi[which(indels_library$psi > 100)] <- 100


for (k in 1:10) {
  y_axis_vector <- c()
  x_axis_start_vector <- c()
  psi_vector <- c()
  
  
  
  for (i in 1:(63-k+1)) {
    
    this_id <- paste("Del_k",
                     k,
                     "_",
                     i,
                     "-",
                     i+k-1,
                     sep = "",
                     collapse = "")
    
    x_axis_start_vector <- c(x_axis_start_vector,
                             i)
    
    if (length(y_axis_vector)==0) {
      y_axis_vector <- c(y_axis_vector,
                         k)
    } else {
      if (y_axis_vector[length(y_axis_vector)]==1) {
        y_axis_vector <- c(y_axis_vector,
                           k)
      } else {
        y_axis_vector <- c(y_axis_vector,
                           y_axis_vector[length(y_axis_vector)]-1)
      }
    }
    
    
    if (any(grepl(pattern = this_id, x = indels_library$id))) {
      psi_vector <- c(psi_vector,
                      indels_library$psi[grep(pattern = this_id, x = indels_library$id)])
    } else {
      psi_vector <- c(psi_vector,
                      NA)
    }
    
  }
  
  
  
  
  plot_tibble <- tibble(y = y_axis_vector,
                        x = x_axis_start_vector + k/2 - 0.5,
                        width = k,
                        psi = psi_vector)
  
  
  
  
  
  colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 7, palette = "Spectral")[2:6]))
  my_colours <- colour_function(11)
  
  ggplot(data = plot_tibble,
         mapping = aes(x=as.factor(as.numeric(x)),
                       y=y,
                       fill=psi,
                       width = width)) +
    geom_tile() +
    scale_fill_gradientn(colours = my_colours, 
                         values = rescale(c(0,
                                            10,
                                            20,
                                            30,
                                            40,
                                            50,
                                            60,
                                            70,
                                            80,
                                            90,
                                            100)),
                         limits = c(0,100),
                         guide = "colourbar",
                         na.value = "black",
                         oob = scales::squish) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          aspect.ratio = k/63,
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    ylab("") +
    xlab("position")
  
  ggsave(filename = paste("132_wallplot_k", k, ".pdf", sep = "", collapse = ""),
         height = 10,
         width = 10,
         useDingbats = FALSE)
}


















for (k in 1:10) {
  print(k)
  psi_vector <- c()
  
  
  
  for (i in 1:(63-k+1)) {
    
    this_id <- paste("Del_k",
                     k,
                     "_",
                     i,
                     "-",
                     i+k-1,
                     sep = "",
                     collapse = "")
    
    
    
    if (any(grepl(pattern = this_id, x = indels_library$id))) {
      psi_vector <- c(psi_vector,
                      indels_library$psi[grep(pattern = this_id, x = indels_library$id)])
    } else {
      psi_vector <- c(psi_vector,
                      NA)
    }
    
  }
  
  
  
  
  print(mean(abs(psi_vector - 49.1) > 10, na.rm = TRUE))
  print(mean(abs(psi_vector - 49.1), na.rm = TRUE))
  print(mean((psi_vector - 49.1) > 10, na.rm = TRUE))
  print(mean((psi_vector - 49.1) < -10, na.rm = TRUE))
  
  
}





















psi_vector <- c()



for (k in 1:10) {
  print(k)
  
  
  
  for (i in 1:(63-k+1)) {
    
    this_id <- paste("Del_k",
                     k,
                     "_",
                     i,
                     "-",
                     i+k-1,
                     sep = "",
                     collapse = "")
    
    
    
    if (any(grepl(pattern = this_id, x = indels_library$id))) {
      psi_vector <- c(psi_vector,
                      indels_library$psi[grep(pattern = this_id, x = indels_library$id)])
    } else {
      psi_vector <- c(psi_vector,
                      NA)
    }
    
  }
  
  
  
  
  
}





print(mean(abs(psi_vector - 49.1) > 10, na.rm = TRUE))
print(mean(abs(psi_vector - 49.1), na.rm = TRUE))
print(mean((psi_vector - 49.1) > 10, na.rm = TRUE))
print(mean((psi_vector - 49.1) < -10, na.rm = TRUE))

