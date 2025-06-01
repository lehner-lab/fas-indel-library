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

indels_library$psi <- (49.1/1) * exp(indels_library$es)
indels_library$psi <- sapply(X = indels_library$psi,
                             FUN = function(x){
                               if (x > 100) {
                                 x <- 100
                               }
                               x
                             })


indels_library <- indels_library %>%
  filter(id != "WT")


write_delim(x = indels_library,
            file = "106_dataset.tsv",
            delim = "\t")


plot_tibble <- tibble(id = indels_library$id,
                      psi = sapply(X = indels_library$psi,
                                   FUN = function(x){
                                     if (x > 100) {
                                       x <- 100
                                     }
                                     x
                                   }),
                      mutation_class = sapply(X = indels_library$sequence,
                                              FUN = function(x) {
                                                result <- "other"
                                                if (nchar(x) < 63) {
                                                  result <- "deletion"
                                                }
                                                
                                                if (nchar(x) == 63) {
                                                  result <- "substitution"
                                                }
                                                
                                                if (nchar(x) > 63) {
                                                  result <- "insertion"
                                                }
                                                
                                                result
                                              }),
                      sequence_length = sapply(X = indels_library$sequence,
                                               FUN = nchar))



plot_tibble <- separate_rows(plot_tibble, id, sep = ";")


ggplot(data = plot_tibble %>% filter(sequence_length > 61 & sequence_length < 65),
       mapping = aes(x = factor(mutation_class,
                                levels = c("substitution","deletion",  "insertion")),
                     y = psi,
                     group = mutation_class)) +
  geom_hline(yintercept = 49.1,
             size = 2,
             colour = "gray60") +
  geom_violin(fill = "gray90", colour = NA, scale = "width") +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  ylab("percent spliced in") +
  xlab("") +
  coord_cartesian(ylim = c(0,100))


write_delim(x = plot_tibble,
            file = "106_figure_01d.tsv",
            delim = "\t")

colnames(plot_tibble)
mean(abs((plot_tibble %>%
            filter(sequence_length == 62))$psi - 49.1) > 10)
mean(abs((plot_tibble %>%
            filter(sequence_length == 62))$psi - 49.1))
mean(((plot_tibble %>%
        filter(sequence_length == 62))$psi - 49.1) > 10)
mean(((plot_tibble %>%
        filter(sequence_length == 62))$psi - 49.1) < -10)


colnames(plot_tibble)
mean(abs((plot_tibble %>%
            filter(mutation_class == "substitution"))$psi - 49.1) > 10)
mean(abs((plot_tibble %>%
            filter(mutation_class == "substitution"))$psi - 49.1))
mean(((plot_tibble %>%
         filter(mutation_class == "substitution"))$psi - 49.1) > 10)
mean(((plot_tibble %>%
         filter(mutation_class == "substitution"))$psi - 49.1) < -10)

colnames(plot_tibble)
mean(abs((plot_tibble %>%
            filter(sequence_length == 64))$psi - 49.1) > 10)
mean(abs((plot_tibble %>%
            filter(sequence_length == 64))$psi - 49.1))
mean(((plot_tibble %>%
         filter(sequence_length == 64))$psi - 49.1) > 10)
mean(((plot_tibble %>%
         filter(sequence_length == 64))$psi - 49.1) < -10)



nrow(plot_tibble %>%
       filter(sequence_length == 64))

nrow(plot_tibble %>%
       filter(sequence_length == 62))

ggsave(filename = "106_distribution_psi_by_mutation_class.pdf",
       height = 10,
       width = 10,
       useDingbats = F)





tibble_1nt <- plot_tibble %>% filter(sequence_length > 61 & sequence_length < 65)

mean(abs(tibble_1nt$psi[which(tibble_1nt$mutation_class == "substitution")] - 49.1) >= 10)
mean(tibble_1nt$psi[which(tibble_1nt$mutation_class == "substitution")])
range(tibble_1nt$psi[which(tibble_1nt$mutation_class == "substitution")])


mean(abs(tibble_1nt$psi[which(tibble_1nt$mutation_class == "deletion")] - 49.1) >= 10)
mean(tibble_1nt$psi[which(tibble_1nt$mutation_class == "deletion")])
range(tibble_1nt$psi[which(tibble_1nt$mutation_class == "deletion")])


mean(abs(tibble_1nt$psi[which(tibble_1nt$mutation_class == "insertion")] - 49.1) >= 10)
mean(tibble_1nt$psi[which(tibble_1nt$mutation_class == "insertion")])
range(tibble_1nt$psi[which(tibble_1nt$mutation_class == "insertion")])







substitutions_vector <- mean(abs(49.1 - (indels_library %>%
                                           filter(nchar(sequence)==63) %>%
                                           filter(id != "WT"))$psi) > 10)


deletions_vector <- c(mean(abs(49.1 - (indels_library %>%
                                         filter(nchar(sequence)<63) %>%
                                         filter(id != "WT"))$psi) > 10),
                      sapply(X = 1:21,
                             FUN = function(x){
                               absolute_delta_psi <- abs(49.1 - (indels_library %>%
                                                                   filter(nchar(sequence) == (63-x) ))$psi)
                               mean(absolute_delta_psi > 10)
                             })
)


insertions_vector <- c(mean(abs(49.1 - (indels_library %>%
                                          filter(nchar(sequence)>63) %>%
                                          filter(id != "WT"))$psi) > 10),
                       sapply(X = 1:4,
                              FUN = function(x){
                                absolute_delta_psi <- abs(49.1 - (indels_library %>%
                                                                    filter(nchar(sequence) == (63+x) ))$psi)
                                mean(absolute_delta_psi > 10)
                              }))

pdf(file = "106_FIGURES_percent_mutations_change_psi_more_than_10units.pdf",
    height = 6,
    width = 10,
    useDingbats = F)
barplot(height = 100*c(substitutions_vector,
                       deletions_vector,
                       insertions_vector),
        border = NA,
        col = c("gray80",
                "gray80",
                rep("gray90", length(deletions_vector)-1),
                "gray80",
                rep("gray90", length(insertions_vector)-1)),
        las = 2,
        names.arg = c("substitutions",
                      "all",
                      as.character(1:20),
                      ">=21",
                      "all",
                      as.character(1:4)), 
        ylab = "percent")
dev.off()










