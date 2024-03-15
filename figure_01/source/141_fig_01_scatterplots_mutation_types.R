
library(tidyverse)
library(scales)
load("data/indels_library.RData")



idx_singles <- which(nchar(indels_library$sequence)==63 & indels_library$id != "WT")
density_object <- density(x = indels_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]

indels_library$es <- indels_library$es - library_bias

indels_library$psi <- (49.1/1) * exp(indels_library$es)

indels_library$psi_min <- (49.1/1) * exp(indels_library$es - 1.96*indels_library$standard.error)

indels_library$psi_max <- (49.1/1) * exp(indels_library$es + 1.96*indels_library$standard.error)


psi_corrections <- indels_library$psi[which(indels_library$psi > 100)] - 100

indels_library$psi_min[which(indels_library$psi > 100)] <- indels_library$psi_min[which(indels_library$psi > 100)] - psi_corrections
indels_library$psi_max[which(indels_library$psi > 100)] <- indels_library$psi_max[which(indels_library$psi > 100)] - psi_corrections
indels_library$psi[which(indels_library$psi > 100)] <- 100


# select substitutions
plot_tibble_1 <- indels_library %>%
  filter(nchar(sequence) == 64) %>%
  separate_rows(id, sep = ";")

plot_tibble_1$position <- sapply(X = plot_tibble_1$id,
                                 FUN = function(x){
                                   as.numeric(strsplit(x = x,
                                                       split = "_")[[1]][3])
                                 })

plot_tibble_1$del_id = sapply(X = plot_tibble_1$position,
                              FUN = function(x){
                                paste("Del_k1_",
                                      x,
                                      "-",
                                      x,
                                      sep = "",
                                      collapse = "")
                              })

plot_tibble_1$inserted_base <- sapply(X = plot_tibble_1$id,
                                      FUN = function(x){
                                        strsplit(x = x, split = "_")[[1]][4]
                                      })


wt_sequence <- indels_library$sequence[indels_library$id == "WT"]

plot_tibble_1$deleted_base <- sapply(X = plot_tibble_1$del_id,
                                     FUN = function(x){
                                       position <- as.numeric(strsplit(x = x, split = "-")[[1]][2])
                                       substr(x = wt_sequence,
                                              start = position,
                                              stop = position)
                                     })


plot_tibble_1$del_es = sapply(X = plot_tibble_1$position,
                              FUN = function(x){
                                
                                this_id <- paste("Del_k1_",
                                                 x,
                                                 "-",
                                                 x,
                                                 sep = "",
                                                 collapse = "")
                                
                                idx <- grep(pattern = this_id,
                                            indels_library$id)
                                
                                if (length(idx) > 0) {
                                  indels_library$es[idx]
                                } else {
                                  NA
                                }
                                
                              })

plot_tibble_1$del_se = sapply(X = plot_tibble_1$position,
                              FUN = function(x){
                                
                                this_id <- paste("Del_k1_",
                                                 x,
                                                 "-",
                                                 x,
                                                 sep = "",
                                                 collapse = "")
                                
                                idx <- grep(pattern = this_id,
                                            indels_library$id)
                                
                                if (length(idx) > 0) {
                                  indels_library$standard.error[idx]
                                } else {
                                  NA
                                }
                                
                                
                              })

plot_tibble_1$del_psi = sapply(X = plot_tibble_1$position,
                               FUN = function(x){
                                 
                                 this_id <- paste("Del_k1_",
                                                  x,
                                                  "-",
                                                  x,
                                                  sep = "",
                                                  collapse = "")
                                 
                                 idx <- grep(pattern = this_id,
                                             indels_library$id)
                                 
                                 if (length(idx) > 0) {
                                   indels_library$psi[idx]
                                 } else {
                                   NA
                                 }
                                 
                               })

plot_tibble_1$del_psi_min = sapply(X = plot_tibble_1$position,
                                   FUN = function(x){
                                     
                                     this_id <- paste("Del_k1_",
                                                      x,
                                                      "-",
                                                      x,
                                                      sep = "",
                                                      collapse = "")
                                     
                                     idx <- grep(pattern = this_id,
                                                 indels_library$id)
                                     
                                     if (length(idx) > 0) {
                                       indels_library$psi_min[idx]
                                     } else {
                                       NA
                                     }
                                     
                                   })

plot_tibble_1$del_psi_max = sapply(X = plot_tibble_1$position,
                                   FUN = function(x){
                                     
                                     this_id <- paste("Del_k1_",
                                                      x,
                                                      "-",
                                                      x,
                                                      sep = "",
                                                      collapse = "")
                                     
                                     idx <- grep(pattern = this_id,
                                                 indels_library$id)
                                     
                                     if (length(idx) > 0) {
                                       indels_library$psi_max[idx]
                                     } else {
                                       NA
                                     }
                                     
                                   })


plot_tibble_1$insertion_location <- "after"








plot_tibble_2 <- indels_library %>%
  filter(nchar(sequence) == 64) %>%
  separate_rows(id, sep = ";")

plot_tibble_2$position <- sapply(X = plot_tibble_2$id,
                                 FUN = function(x){
                                   as.numeric(strsplit(x = x,
                                                       split = "_")[[1]][3])
                                 })
plot_tibble_2$position <- plot_tibble_2$position - 1

plot_tibble_2$del_id = sapply(X = plot_tibble_2$position,
                              FUN = function(x){
                                paste("Del_k1_",
                                      x,
                                      "-",
                                      x,
                                      sep = "",
                                      collapse = "")
                              })

plot_tibble_2$inserted_base <- sapply(X = plot_tibble_2$id,
                                      FUN = function(x){
                                        strsplit(x = x, split = "_")[[1]][4]
                                      })


wt_sequence <- indels_library$sequence[indels_library$id == "WT"]

plot_tibble_2$deleted_base <- sapply(X = plot_tibble_2$del_id,
                                     FUN = function(x){
                                       position <- as.numeric(strsplit(x = x, split = "-")[[1]][2])
                                       substr(x = wt_sequence,
                                              start = position,
                                              stop = position)
                                     })


plot_tibble_2$del_es = sapply(X = plot_tibble_2$position,
                              FUN = function(x){
                                
                                this_id <- paste("Del_k1_",
                                                 x,
                                                 "-",
                                                 x,
                                                 sep = "",
                                                 collapse = "")
                                
                                idx <- grep(pattern = this_id,
                                            indels_library$id)
                                
                                if (length(idx) > 0) {
                                  indels_library$es[idx]
                                } else {
                                  NA
                                }
                                
                              })

plot_tibble_2$del_se = sapply(X = plot_tibble_2$position,
                              FUN = function(x){
                                
                                this_id <- paste("Del_k1_",
                                                 x,
                                                 "-",
                                                 x,
                                                 sep = "",
                                                 collapse = "")
                                
                                idx <- grep(pattern = this_id,
                                            indels_library$id)
                                
                                if (length(idx) > 0) {
                                  indels_library$standard.error[idx]
                                } else {
                                  NA
                                }
                                
                                
                              })

plot_tibble_2$del_psi = sapply(X = plot_tibble_2$position,
                               FUN = function(x){
                                 
                                 this_id <- paste("Del_k1_",
                                                  x,
                                                  "-",
                                                  x,
                                                  sep = "",
                                                  collapse = "")
                                 
                                 idx <- grep(pattern = this_id,
                                             indels_library$id)
                                 
                                 if (length(idx) > 0) {
                                   indels_library$psi[idx]
                                 } else {
                                   NA
                                 }
                                 
                               })

plot_tibble_2$del_psi_min = sapply(X = plot_tibble_2$position,
                                   FUN = function(x){
                                     
                                     this_id <- paste("Del_k1_",
                                                      x,
                                                      "-",
                                                      x,
                                                      sep = "",
                                                      collapse = "")
                                     
                                     idx <- grep(pattern = this_id,
                                                 indels_library$id)
                                     
                                     if (length(idx) > 0) {
                                       indels_library$psi_min[idx]
                                     } else {
                                       NA
                                     }
                                     
                                   })

plot_tibble_2$del_psi_max = sapply(X = plot_tibble_2$position,
                                   FUN = function(x){
                                     
                                     this_id <- paste("Del_k1_",
                                                      x,
                                                      "-",
                                                      x,
                                                      sep = "",
                                                      collapse = "")
                                     
                                     idx <- grep(pattern = this_id,
                                                 indels_library$id)
                                     
                                     if (length(idx) > 0) {
                                       indels_library$psi_max[idx]
                                     } else {
                                       NA
                                     }
                                     
                                   })

plot_tibble_2$insertion_location <- "before"


plot_tibble_insertions_deletions <- rbind(plot_tibble_1,
                                          plot_tibble_2)






single_mutants <- indels_library$id[idx_singles]


plot_tibble_subs <- tibble(position = sapply(X = single_mutants,
                                             FUN = function(x){
                                               as.numeric(substr(x = x,
                                                                 start = 2,
                                                                 stop = nchar(x)-1))
                                             }),
                           
                           sub_id = single_mutants,
                           
                           sub_es = sapply(X = single_mutants,
                                           FUN = function(x){
                                             
                                             idx <- grep(pattern = x,
                                                         indels_library$id)
                                             
                                             if (length(idx) > 0) {
                                               indels_library$es[idx]
                                             } else {
                                               NA
                                             }
                                             
                                           }),
                           
                           sub_se = sapply(X = single_mutants,
                                           FUN = function(x){
                                             
                                             idx <- grep(pattern = x,
                                                         indels_library$id)
                                             
                                             if (length(idx) > 0) {
                                               indels_library$standard.error[idx]
                                             } else {
                                               NA
                                             }
                                             
                                             
                                           }),
                           
                           sub_psi = sapply(X = single_mutants,
                                            FUN = function(x){
                                              
                                              idx <- grep(pattern = x,
                                                          indels_library$id)
                                              
                                              if (length(idx) > 0) {
                                                indels_library$psi[idx]
                                              } else {
                                                NA
                                              }
                                              
                                            }),
                           
                           sub_psi_min = sapply(X = single_mutants,
                                                FUN = function(x){
                                                  
                                                  idx <- grep(pattern = x,
                                                              indels_library$id)
                                                  
                                                  if (length(idx) > 0) {
                                                    indels_library$psi_min[idx]
                                                  } else {
                                                    NA
                                                  }
                                                  
                                                }),
                           
                           sub_psi_max = sapply(X = single_mutants,
                                                FUN = function(x){
                                                  
                                                  idx <- grep(pattern = x,
                                                              indels_library$id)
                                                  
                                                  if (length(idx) > 0) {
                                                    indels_library$psi_max[idx]
                                                  } else {
                                                    NA
                                                  }
                                                  
                                                }))




plot_tibble <- merge(x = plot_tibble_insertions_deletions,
                     y = plot_tibble_subs,
                     by = "position")



?distinct

# substitutions vs deletions


my_correlation <- signif(x = as.numeric(cor.test(x = (plot_tibble %>%
                                                        select(sub_psi, del_psi) %>%
                                                        distinct())$sub_psi,
                                                 y = (plot_tibble %>%
                                                        select(sub_psi, del_psi) %>%
                                                        distinct())$del_psi,
                                                 method = "spearman")$estimate),
                         digits = 2)


ggplot(data = plot_tibble %>%
         select(sub_psi, sub_psi_max, sub_psi_min,
                del_psi, del_psi_max, del_psi_min,
                position) %>%
         distinct(),
       mapping = aes(x = sub_psi,
                     y = del_psi)) +
  # geom_abline(slope = 1,
  #             intercept = 0) +
  geom_vline(xintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_hline(yintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_errorbar(aes(ymin=del_psi_min,
                    ymax=del_psi_max),
                colour = "gray30") +
  geom_errorbar(aes(xmin = sub_psi_min,
                    xmax = sub_psi_max),
                colour = "gray30") +
  geom_point(size = 3,
             colour = "gray30") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  ylab("PSI (deletion)") +
  xlab("PSI (substitution)") +
  coord_cartesian(xlim = c(0,90),
                  ylim = c(0,90)) +
  ggtitle(paste("rho = ",
                my_correlation,
                sep = "",
                collapse = ""))


ggsave(filename = "141_deletion_vs_substitutions.pdf",
       width = 7,
       height = 7,
       useDingbats = F)



# deletions vs insertions


my_correlation <- signif(x = as.numeric(cor.test(x = (plot_tibble %>%
                                                        select(psi, psi_max, psi_min,
                                                               del_psi, del_psi_max, del_psi_min,
                                                               insertion_location) %>%
                                                        filter(insertion_location == "after") %>%
                                                        distinct())$psi,
                                                 y = (plot_tibble %>%
                                                        select(psi, psi_max, psi_min,
                                                               del_psi, del_psi_max, del_psi_min,
                                                               insertion_location) %>%
                                                        filter(insertion_location == "after") %>%
                                                        distinct())$del_psi,
                                                 method = "spearman")$estimate),
                         digits = 2)


ggplot(data = plot_tibble %>%
         select(psi, psi_max, psi_min,
                del_psi, del_psi_max, del_psi_min,
                insertion_location) %>%
         filter(insertion_location == "after") %>%
         distinct(),
       mapping = aes(x = psi,
                     y = del_psi)) +
  # geom_abline(slope = 1,
  #             intercept = 0) +
  geom_vline(xintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_hline(yintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_errorbar(aes(ymin=del_psi_min,
                    ymax=del_psi_max),
                colour = "gray30") +
  geom_errorbar(aes(xmin = psi_min,
                    xmax = psi_max),
                colour = "gray30") +
  geom_point(size = 3,
             colour = "gray30") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  ylab("PSI (insertion after)") +
  xlab("PSI (deletion)") +
  coord_cartesian(xlim = c(0,100),
                  ylim = c(0,100)) +
  ggtitle(paste("rho = ",
                my_correlation,
                sep = "",
                collapse = ""))

ggsave(filename = "141_deletion_vs_insertionsAfter.pdf",
       width = 7,
       height = 7,
       useDingbats = F)





my_correlation <- signif(x = as.numeric(cor.test(x = (plot_tibble %>%
                                                        select(psi, psi_max, psi_min,
                                                               del_psi, del_psi_max, del_psi_min,
                                                               insertion_location) %>%
                                                        filter(insertion_location == "before") %>%
                                                        distinct())$psi,
                                                 y = (plot_tibble %>%
                                                        select(psi, psi_max, psi_min,
                                                               del_psi, del_psi_max, del_psi_min,
                                                               insertion_location) %>%
                                                        filter(insertion_location == "before") %>%
                                                        distinct())$del_psi,
                                                 method = "spearman")$estimate),
                         digits = 2)


ggplot(data = plot_tibble %>%
         select(psi, psi_max, psi_min,
                del_psi, del_psi_max, del_psi_min,
                insertion_location) %>%
         filter(insertion_location == "before") %>%
         distinct(),
       mapping = aes(x = psi,
                     y = del_psi)) +
  # geom_abline(slope = 1,
  #             intercept = 0) +
  geom_vline(xintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_hline(yintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_errorbar(aes(ymin=del_psi_min,
                    ymax=del_psi_max),
                colour = "gray30") +
  geom_errorbar(aes(xmin = psi_min,
                    xmax = psi_max),
                colour = "gray30") +
  geom_point(size = 3,
             colour = "gray30") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  ylab("PSI (insertion before)") +
  xlab("PSI (deletion)") +
  coord_cartesian(xlim = c(0,100),
                  ylim = c(0,100)) +
  ggtitle(paste("rho = ",
                my_correlation,
                sep = "",
                collapse = ""))

ggsave(filename = "141_deletion_vs_insertionsBefore.pdf",
       width = 7,
       height = 7,
       useDingbats = F)







# substitutions vs insertions

my_correlation <- signif(x = as.numeric(cor.test(x = (plot_tibble %>%
                                                        select(psi, psi_max, psi_min,
                                                               sub_psi, sub_psi_max, sub_psi_min,
                                                               insertion_location) %>%
                                                        filter(insertion_location == "after") %>%
                                                        distinct())$psi,
                                                 y = (plot_tibble %>%
                                                        select(psi, psi_max, psi_min,
                                                               sub_psi, sub_psi_max, sub_psi_min,
                                                               insertion_location) %>%
                                                        filter(insertion_location == "after") %>%
                                                        distinct())$sub_psi,
                                                 method = "spearman")$estimate),
                         digits = 2)

ggplot(data = plot_tibble %>%
         select(psi, psi_max, psi_min,
                sub_psi, sub_psi_max, sub_psi_min,
                insertion_location) %>%
         filter(insertion_location == "after") %>%
         distinct(),
       mapping = aes(x = psi,
                     y = sub_psi)) +
  # geom_abline(slope = 1,
  #             intercept = 0) +
  geom_vline(xintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_hline(yintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_errorbar(aes(ymin=sub_psi_min,
                    ymax=sub_psi_max),
                colour = "gray30") +
  geom_errorbar(aes(xmin = psi_min,
                    xmax = psi_max),
                colour = "gray30") +
  geom_point(size = 3,
             colour = "gray30") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  ylab("PSI (insertion after)") +
  xlab("PSI (substitution)") +
  coord_cartesian(xlim = c(0,100),
                  ylim = c(0,100)) +
  ggtitle(paste("rho = ",
                my_correlation,
                sep = "",
                collapse = ""))

ggsave(filename = "141_substitutions_vs_insertionsAfter.pdf",
       width = 7,
       height = 7,
       useDingbats = F)





my_correlation <- signif(x = as.numeric(cor.test(x = (plot_tibble %>%
                                                        select(psi, psi_max, psi_min,
                                                               sub_psi, sub_psi_max, sub_psi_min,
                                                               insertion_location) %>%
                                                        filter(insertion_location == "before") %>%
                                                        distinct())$psi,
                                                 y = (plot_tibble %>%
                                                        select(psi, psi_max, psi_min,
                                                               sub_psi, sub_psi_max, sub_psi_min,
                                                               insertion_location) %>%
                                                        filter(insertion_location == "before") %>%
                                                        distinct())$sub_psi,
                                                 method = "spearman")$estimate),
                         digits = 2)

ggplot(data = plot_tibble %>%
         select(psi, psi_max, psi_min,
                sub_psi, sub_psi_max, sub_psi_min,
                insertion_location) %>%
         filter(insertion_location == "before") %>%
         distinct(),
       mapping = aes(x = psi,
                     y = sub_psi)) +
  # geom_abline(slope = 1,
  #             intercept = 0) +
  geom_vline(xintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_hline(yintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_errorbar(aes(ymin=sub_psi_min,
                    ymax=sub_psi_max),
                colour = "gray30") +
  geom_errorbar(aes(xmin = psi_min,
                    xmax = psi_max),
                colour = "gray30") +
  geom_point(size = 3,
             colour = "gray30") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  ylab("PSI (insertion before)") +
  xlab("PSI (substitution)") +
  coord_cartesian(xlim = c(0,100),
                  ylim = c(0,100)) +
  ggtitle(paste("rho = ",
                my_correlation,
                sep = "",
                collapse = ""))

ggsave(filename = "141_substitutions_vs_insertionsBefore.pdf",
       width = 7,
       height = 7,
       useDingbats = F)








insertions_tibble <- merge(x = (plot_tibble %>%
                                  filter(insertion_location == "after") %>%
                                  select(position, psi, psi_max, psi_min, id) %>%
                                  distinct()),
                           y = (plot_tibble %>%
                                  filter(insertion_location == "before") %>%
                                  select(position, psi, psi_max, psi_min, id) %>%
                                  distinct()),
                           by = "position")


my_correlation <- signif(x = as.numeric(cor.test(x = insertions_tibble$psi.x,
                                                 y = insertions_tibble$psi.y,
                                                 method = "spearman")$estimate),
                         digits = 2)

ggplot(data = insertions_tibble,
       mapping = aes(x = psi.x,
                     y = psi.y)) +
  # geom_abline(slope = 1,
  #             intercept = 0) +
  geom_vline(xintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_hline(yintercept = 49.1,
             lty = 2,
             colour = "gray85") +
  geom_errorbar(aes(ymin = psi_min.y,
                    ymax = psi_max.y),
                colour = "gray30") +
  geom_errorbar(aes(xmin = psi_min.x,
                    xmax = psi_max.x),
                colour = "gray30") +
  geom_point(size = 3,
             colour = "gray30") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  ylab("PSI (insertion before)") +
  xlab("PSI (insertion after)") +
  coord_cartesian(xlim = c(0,100),
                  ylim = c(0,100)) +
  ggtitle(paste("rho = ",
                my_correlation,
                sep = "",
                collapse = ""))

ggsave(filename = "141_insertionsAfter_vs_insertionsBefore.pdf",
       width = 7,
       height = 7,
       useDingbats = F)



