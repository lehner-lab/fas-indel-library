
library(tidyverse)
library(scales)
load("data/indels_library.RData")
load("model_predictions/spliceAI/spliceAI_predictions.RData")

wt_sequence <- indels_library$sequence[which(indels_library$id == "WT")]
idx_singles <- which(nchar(indels_library$sequence)==63 & indels_library$id != "WT")
density_object <- density(x = indels_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]

indels_library$es <- indels_library$es - library_bias

indels_library$psi <- (49.1/1) * exp(indels_library$es)

indels_library$psi[which(indels_library$psi>100)] <-  100


indels_library$spliceai <- sapply(X = indels_library$id,
                                  FUN = function(x){
                                    id <- strsplit(x, ";")[[1]][1]
                                    if (id == "WT") {
                                      NA
                                    } else {
                                      results_table$spliceAI[which(results_table$id == id)]
                                    }
                                  })



psi_vector <- c()
spliceAI_vector <- c()
position_vector <- c()
k_vector <- c()

# sliding windows size 1-10
# each window deletion is defined by the one and only deletion that exists there
for (k in 1:10) {
  
  # each position
  for (p in 1:(63-k+1)) {
    
    this_id <- paste("Del_k",
                     k,
                     "_",
                     p,
                     "-",
                     p+k-1,
                     sep = "",
                     collapse = "")
    
    idx <- grep(pattern = this_id, x = indels_library$id)
    
    if (length(idx)>0) {
      psi_vector <- c(psi_vector,
                      indels_library$psi[idx])
      spliceAI_vector <- c(spliceAI_vector,
                           indels_library$spliceai[idx])
      position_vector <- c(position_vector,
                           p + (k/2) - 0.5)
      k_vector <- c(k_vector,
                    k)
    }
    
    
  }
  
  
}





plot_tibble <- tibble(psi = psi_vector,
                      spliceAI = spliceAI_vector,
                      k = k_vector,
                      position = position_vector)


ggplot(data = plot_tibble,
       mapping = aes(x = position,
                     y = psi,
                     colour = as.factor(k))) +
  geom_hline(yintercept = 49.1,
             colour = "gray80",
             size = 2) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank())



colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 11, palette = "Spectral")[2:10]))
my_colours <- colour_function(11)


ggplot(data = plot_tibble %>% filter(k %in% 1:6),
       mapping = aes(x = position,
                     y = spliceAI)) +
  geom_hline(yintercept = 0,
             colour = "gray90",
             size = 2) +
  geom_point(colour = "gray30",
             fill = "gray30",
             mapping = aes(shape = as.factor(k))) +
  geom_smooth(method = "loess",
              span = 0.15,
              colour = my_colours[3],
              fill = my_colours[3],
              alpha = 0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 20/63) +
  scale_x_continuous(limits = c(1,63),
                     expand = c(0.01, 0.01),
                     breaks = 1:63) +
  scale_shape_manual(values = c(21:25,3),
                     name = "deletion length")


ggsave(filename = "136B_exon_modules_plot.pdf",
       width = 10,
       height = 10,
       useDingbats = F)



