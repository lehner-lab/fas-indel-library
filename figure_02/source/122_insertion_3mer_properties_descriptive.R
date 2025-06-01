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

# 1mer insertions


k <- 3
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



# rowclust <- hclust(dist(scores_matrix))
# ordered_kmers <- rownames(scores_matrix)[rowclust$order]



list_a <- vector(mode = "list", length = 4)
list_c <- vector(mode = "list", length = 4)
list_g <- vector(mode = "list", length = 4)
list_t <- vector(mode = "list", length = 4)

for (each_3mer in rownames(scores_matrix)) {
  counts_a <- str_count(string = each_3mer, pattern = "A")
  counts_c <- str_count(string = each_3mer, pattern = "C")
  counts_g <- str_count(string = each_3mer, pattern = "G")
  counts_t <- str_count(string = each_3mer, pattern = "T")
  
  list_a[[counts_a+1]] <- c(list_a[[counts_a+1]],
                            scores_matrix[each_3mer,])
  list_c[[counts_c+1]] <- c(list_c[[counts_c+1]],
                            scores_matrix[each_3mer,])
  list_g[[counts_g+1]] <- c(list_g[[counts_g+1]],
                            scores_matrix[each_3mer,])
  list_t[[counts_t+1]] <- c(list_t[[counts_t+1]],
                            scores_matrix[each_3mer,])
}

names(list_a) <- c(0, 1, 2, 3)
names(list_c) <- c(0, 1, 2, 3)
names(list_g) <- c(0, 1, 2, 3)
names(list_t) <- c(0, 1, 2, 3)



plot_tibble <- tibble(nucleotide = c(rep("A",
                                         length(unlist(list_a))),
                                     rep("C",
                                         length(unlist(list_c))),
                                     rep("G",
                                         length(unlist(list_g))),
                                     rep("T",
                                         length(unlist(list_t)))),
                      counts = c(rep(0:3, sapply(list_a, length)),
                                 rep(0:3, sapply(list_c, length)),
                                 rep(0:3, sapply(list_g, length)),
                                 rep(0:3, sapply(list_t, length))),
                      psi = c(unlist(list_a),
                             unlist(list_c),
                             unlist(list_g),
                             unlist(list_t))
                      )

write_delim(x = plot_tibble,
            delim = "\t",
            file = "122_figure02e.tsv")

for (each_nucleotide in c("A", "C", "G", "T")) {
  ggplot(data = plot_tibble %>%
           filter(nucleotide == each_nucleotide),
         mapping = aes(x = counts,
                       y = psi,
                       group = counts)) +
    geom_violin(scale = "width",
                trim = T,
                fill = "gray90") +
    geom_boxplot(width = 0.1,
                 fill = "white",
                 outlier.shape = NA) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          aspect.ratio = 1) +
    xlab("counts") +
    ylab("percent spliced in") +
    coord_cartesian(ylim = c(0,100)) +
    ggtitle(each_nucleotide)
  
  file_name <- paste("122_3mers_",
                     each_nucleotide,
                     "_content.pdf",
                     sep = "",
                     collapse = "")
  
  ggsave(filename = file_name,
         height = 10,
         width = 10,
         useDingbats = FALSE)
  
  
  
}





