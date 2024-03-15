
library(tidyverse)
library(scales)
load("data/indels_library.RData")



idx_singles <- which(nchar(indels_library$sequence)==63 & indels_library$id != "WT")
density_object <- density(x = indels_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]

indels_library$es <- indels_library$es - library_bias

indels_library$psi <- (49.1/1) * exp(indels_library$es)


# 1mer insertions


k <- 2
idx_k <- which(nchar(indels_library$sequence)==63+k)

possible_insertions <- sort(apply(X = expand.grid(rep(list(c("A","T","G","C")), k)),
                                  MARGIN = 1,
                                  FUN = function(x){
                                    paste(x, sep = "", collapse = "")
                                  }))


scores_matrix <- matrix(data = NA,
                        nrow = length(possible_insertions),
                        ncol = 62,
                        dimnames = list(possible_insertions,
                                        1:62))

all_possible_ids <- apply(X = expand.grid(1:62,
                                          possible_insertions),
                          MARGIN = 1,
                          FUN = function(x){
                            paste("Ins_k",
                                  k,
                                  "_",
                                  as.numeric(x[1]),
                                  "_",
                                  x[2],
                                  sep = "",
                                  collapse = "")
                          })


for (each_id in all_possible_ids) {
  idx_id <- grep(pattern = each_id,
                 x = indels_library$id)
  
  if (length(idx_id)==0){
    next
  }
  
  position <- as.numeric(strsplit(x = each_id,
                                  split = "_")[[1]][3])
  insertion <- as.character(strsplit(x = each_id,
                                     split = "_")[[1]][4])
  
  scores_matrix[insertion, position] <- indels_library$psi[idx_id]
  
  
}


scores_matrix[which(scores_matrix>100)] <- 100

order_kmers <- order(apply(X = scores_matrix,
                           MARGIN = 1,
                           FUN = median, na.rm = T))

median_scores <- apply(X = scores_matrix,
                       MARGIN = 1,
                       FUN = median, na.rm = T)
median_scores <- median_scores[rev(order_kmers)]



median_scores_library <- median_scores
save(median_scores_library, file = "110_median_scores.RData")


colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 11, palette = "Spectral")[2:10]))


map2color2<-function(x,pal,breaks){
  pal[findInterval(x,breaks, all.inside=TRUE)]
}
my_palette <- colour_function(100)
my_breaks <- seq(0,100,1)
my_colours <- map2color2(x = median_scores,
                         pal = my_palette,
                         breaks = my_breaks)



scores_matrix %>%
  as.data.frame() %>%
  rownames_to_column("insertion") %>%
  pivot_longer(-c(insertion),
               names_to = "position",
               values_to = "psi") %>%
  ggplot(aes(x=factor(insertion, levels = rownames(scores_matrix)[rev(order_kmers)]),
             y=psi,
             group=insertion)) +
  geom_hline(yintercept = 49.1, colour = "gray80") +
  geom_boxplot(fill = my_colours,
               outlier.size = 0.25) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 35/nrow(scores_matrix)) +
  coord_cartesian(ylim = c(0,100)) +
  xlab("") +
  ylab("percent spliced in")




rownames(scores_matrix) <- gsub(pattern = "T", replacement = "U", x = rownames(scores_matrix))

scores_matrix %>%
  as.data.frame() %>%
  rownames_to_column("insertion") %>%
  pivot_longer(-c(insertion),
               names_to = "position",
               values_to = "psi") %>%
  ggplot(aes(x=factor(insertion, levels = rownames(scores_matrix)[rev(order_kmers)]),
             y=psi,
             group=insertion)) +
  # geom_hline(yintercept = 49.1, colour = "gray80") +
  geom_violin(scale = "width",
              fill = "gray80",
              colour = NA) +
  geom_boxplot(width = 0.4,
               fill = "white",
               outlier.shape = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 35/nrow(scores_matrix)) +
  coord_cartesian(ylim = c(0,100)) +
  xlab("") +
  ylab("percent spliced in")



ggsave(filename = "110_2mer_insertion_boxplots.pdf",
       height = 10,
       width = 10,
       useDingbats = FALSE)



