library(tidyverse)
library(scales)

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


vector_starts <- c()
vector_ends <- c()
vector_es <- c()

for (each_start in 1:63) {
  for (each_end in each_start:63) {
    vector_starts <- c(vector_starts,
                       each_start)
    vector_ends <- c(vector_ends,
                     each_end)
    
    this_k <- each_end - each_start + 1
    
    this_id <- paste("Del_k",
                     this_k,
                     "_",
                     each_start,
                     "-",
                     each_end,
                     sep = "",
                     collapse = "")
    
    idx_this_id <- grep(pattern = this_id,
                        x = indels_library$id)
    
    if (length(idx_this_id)==1) {
      
      vector_es <- c(vector_es,
                     indels_library$es[idx_this_id])
    } else {
      vector_es <- c(vector_es,
                     NA)
    }
    
  }
}



vector_psi <- (49.1/1) * exp(vector_es)

  
colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 7, palette = "Spectral")[2:6]))
my_colours <- colour_function(11)


ggplot(data = tibble(start = vector_starts,
                     end = vector_ends,
                     psi = vector_psi),
       mapping = aes(x = end,
                     y = start,
                     fill = psi)) +
  geom_tile() +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
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
                       oob = scales::squish)

ggsave(filename = "133_figure_02_deletions.pdf",
       height = 10,
       width = 10,
       useDingbats = FALSE)




write_delim(x = tibble(start = vector_starts,
                       end = vector_ends,
                       psi = vector_psi),
            file = "133_figure_05b_dms.tsv",
            delim = "\t")
