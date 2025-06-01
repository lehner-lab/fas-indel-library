setwd("~/CRG/Indel")
library(tidyverse)
library(scales)
load("data/indels_library.RData")

wt_sequence <- indels_library$sequence[which(indels_library$id == "WT")]



spliceAI_table <- read_delim(file = "model_predictions/spliceAI/combined_vcfs.txt", delim = "\t")
spliceAI_table$predictions <- sapply(X = spliceAI_table$splicing,
                                     FUN = function(x){
                                       
                                       delta_scores <- strsplit(x = x,
                                                                split = "|",
                                                                fixed = TRUE)[[1]][3:6]
                                       delta_scores <- as.numeric(delta_scores)
                                       
                                       idx <- which(delta_scores == max(delta_scores))
                                       
                                       result <- delta_scores[idx[1]]
                                       
                                       if (idx[1] %in% c(2,4)) {
                                         result <- result * -1
                                       }
                                       
                                       result
                                       
                                     })


range(spliceAI_table$predictions[grep(pattern = "Del", x = spliceAI_table$id)])




vector_starts <- c()
vector_ends <- c()
vector_delta_psi <- c()

for (each_start in 1:nchar(wt_sequence)) {
  for (each_end in each_start:nchar(wt_sequence)) {
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
                        x = spliceAI_table$id)
    
    if (length(idx_this_id)==1) {
      vector_delta_psi <- c(vector_delta_psi,
                            spliceAI_table$predictions[idx_this_id])
    } else {
      vector_delta_psi <- c(vector_delta_psi,
                            NA)
    }
    
  }
}





colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 7, palette = "Spectral")[2:6]))
my_colours <- colour_function(11)



my_plot <- ggplot(data = tibble(start = vector_starts,
                                end = vector_ends,
                                dPSI = vector_delta_psi),
                  mapping = aes(x = end,
                                y = start,
                                fill = dPSI)) +
  geom_tile() +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  scale_fill_gradientn(colours = my_colours,
                       # values = scales::rescale(seq(-0.25,0.25,length.out = 11)),
                       values = rescale(c(seq(-0.55,-0.1, length.out = 5),
                                          0,
                                          seq(0.05, 0.3, length.out = 5))),
                       limits = c(-0.55,0.3),
                       guide = "colourbar",
                       na.value = "black",
                       oob = scales::squish) +
  coord_flip() +
  scale_x_reverse() +
  ggtitle("spliceAI")



pdf(file = "137B_spliceAI_predictions_mireia_plot.pdf",
    height = 10,
    width = 10)
my_plot
dev.off()


write_delim(x = tibble(start = vector_starts,
                       end = vector_ends,
                       dPSI = vector_delta_psi),
            file = "137B_figure_05b_spliceai.tsv",
            delim = "\t")


