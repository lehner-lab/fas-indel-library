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


k <- 2
idx_k <- which(nchar(indels_library$sequence)==63+k)

possible_insertions <- sort(apply(X = expand.grid(rep(list(c("A","C","G","T")), k)),
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




mean(abs(scores_matrix - 49.1) > 10, na.rm = TRUE)
mean(abs(scores_matrix - 49.1), na.rm = TRUE)
mean((scores_matrix - 49.1) > 10, na.rm = TRUE)
mean((scores_matrix - 49.1) < -10, na.rm = TRUE)







colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 7, palette = "Spectral")[2:6]))
my_colours <- colour_function(11)


scores_matrix %>%
  as.data.frame() %>%
  rownames_to_column("insertion") %>%
  pivot_longer(-c(insertion),
               names_to = "position",
               values_to = "psi") %>%
  ggplot(aes(x=as.factor(as.numeric(position)),
             y=insertion,
             fill=psi)) +
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
        aspect.ratio = length(possible_insertions)/62) +
  scale_y_discrete(limits=rev) +
  ylab("insertion") +
  xlab("position")

ggsave(filename = "135_figure_02_insertions_k2.pdf",
       height = 10,
       width = 10,
       useDingbats = FALSE)













rowclust <- hclust(dist(scores_matrix), method = "ward.D2")
ordered_kmers <- rownames(scores_matrix)[rowclust$order]

pdf(file = "135_insertions_k2_dendrogram.pdf",
    height = 10,
    width = 10,
    useDingbats = FALSE)
plot(rowclust, hang = -1)
dev.off()


scores_matrix %>% 
  as.data.frame() %>%
  rownames_to_column("insertion") %>%
  pivot_longer(-c(insertion),
               names_to = "position",
               values_to = "psi") %>%
  ggplot(aes(x=as.factor(as.numeric(position)),
             y=factor(insertion, levels = ordered_kmers),
             fill=psi)) + 
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
        aspect.ratio = length(possible_insertions)/62) +
  ylab("insertion") +
  xlab("position")

ggsave(filename = "135_figure_02_insertions_k2_CLUSTERED.pdf",
       height = 12,
       width = 10,
       useDingbats = FALSE)
