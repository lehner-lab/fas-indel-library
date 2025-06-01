
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


indels_library$spliceAI <- sapply(X = indels_library$id,
                                  FUN = function(x){
                                    id <- strsplit(x, ";")[[1]][1]
                                    if (id == "WT") {
                                      NA
                                    } else {
                                      results_table$spliceAI[which(results_table$id == id)]
                                    }
                                  })

indels_library$psi_max <- (49.1/1) * exp(indels_library$es + indels_library$standard.error)
indels_library$psi_min <- (49.1/1) * exp(indels_library$es - indels_library$standard.error)





antisense_oligos <- tibble(id = c("FAON_1_-5_16",
                                  "FAON_2_0_21",
                                  "FAON_3_6_26",
                                  "FAON_4_11_31",
                                  "FAON_5_16_36",
                                  "FAON_6_21_41",
                                  "FAON_7_26_46",
                                  "FAON_8_31_51",
                                  "FAON_9_36_56",
                                  "FAON_10_41_61",
                                  "FAON_11_46_+3",
                                  "FAON_12_51_+8"),
                           outside_exon = c(T,
                                            F,F,F,F,F,F,F,F,F,
                                            T,T),
                           corresponding_deletion = c("Del_k16_1-16",
                                                      "Del_k21_1-21",
                                                      "Del_k21_6-26",
                                                      "Del_k21_11-31",
                                                      "Del_k21_16-36",
                                                      "Del_k21_21-41",
                                                      "Del_k21_26-46",
                                                      "Del_k21_31-51",
                                                      "Del_k21_36-56",
                                                      "Del_k21_41-61",
                                                      "Del_k18_46-63",
                                                      "Del_k13_51-63"),
                           psi = c(34.9,
                                   22.0,
                                   33.0,
                                   42.2,
                                   27.4,
                                   41.4,
                                   32.6,
                                   33,
                                   29.1,
                                   33.2,
                                   37.0,
                                   49.3),
                           psi_sd = c(2.3,
                                      1.7,
                                      0.6,
                                      1.3,
                                      0.4,
                                      1.5,
                                      1.6,
                                      1.2,
                                      1.9,
                                      0.8,
                                      2.7,
                                      1.4))











antisense_oligos$del_es <- sapply(X = antisense_oligos$corresponding_deletion,
                                  FUN = function(x){
                                    indels_library$es[grep(pattern = x,
                                                           x = indels_library$id)]
                                  })

antisense_oligos$del_es_se <- sapply(X = antisense_oligos$corresponding_deletion,
                                     FUN = function(x){
                                       indels_library$standard.error[grep(pattern = x,
                                                                          x = indels_library$id)]
                                     })






antisense_oligos$del_psi <- sapply(X = antisense_oligos$corresponding_deletion,
                                   FUN = function(x){
                                     indels_library$psi[grep(pattern = x,
                                                             x = indels_library$id)]
                                   })

antisense_oligos$min_del_psi <- sapply(X = antisense_oligos$corresponding_deletion,
                                       FUN = function(x){
                                         indels_library$psi_min[grep(pattern = x,
                                                                     x = indels_library$id)]
                                       })

antisense_oligos$max_del_psi <- sapply(X = antisense_oligos$corresponding_deletion,
                                       FUN = function(x){
                                         indels_library$psi_max[grep(pattern = x,
                                                                     x = indels_library$id)]
                                       })



antisense_oligos$del_spliceAI <- sapply(X = antisense_oligos$corresponding_deletion,
                                        FUN = function(x){
                                          indels_library$spliceAI[grep(pattern = x,
                                                                       x = indels_library$id)]
                                        })

plot(antisense_oligos$psi,
     exp(antisense_oligos$del_es))

plot(antisense_oligos$psi,
     antisense_oligos$del_spliceAI)
plot(antisense_oligos$psi[which(antisense_oligos$outside_exon == F)],
     (antisense_oligos$del_spliceAI[which(antisense_oligos$outside_exon == F)]))
plot(antisense_oligos$psi[which(antisense_oligos$outside_exon == F)],
     (antisense_oligos$del_psi[which(antisense_oligos$outside_exon == F)]))
abline(0,1)






cor.test(x = antisense_oligos$psi[which(antisense_oligos$outside_exon == F)],
         y = antisense_oligos$del_es[which(antisense_oligos$outside_exon == F)],
         method = "spearman")

cor.test(x = antisense_oligos$psi[which(antisense_oligos$outside_exon == F)],
         y = antisense_oligos$del_spliceAI[which(antisense_oligos$outside_exon == F)],
         method = "spearman")

cor.test(x = antisense_oligos$del_psi[which(antisense_oligos$outside_exon == F)],
         y = antisense_oligos$del_spliceAI[which(antisense_oligos$outside_exon == F)],
         method = "spearman")
plot(x = antisense_oligos$del_psi[which(antisense_oligos$outside_exon == F)],
         y = antisense_oligos$del_spliceAI[which(antisense_oligos$outside_exon == F)])

# 
# ggplot(data = antisense_oligos %>%
#          filter(! outside_exon),
#        mapping = aes(x = psi,
#                      y = exp(del_es))) +
#   geom_smooth(data = subset(antisense_oligos,
#                             outside_exon == F),
#               method = "lm",
#               colour = "black") +
#   geom_point(size = 2) +
#   geom_errorbar(aes(ymin = exp(del_es - 1.96*del_es_se),
#                     ymax = exp(del_es + 1.96*del_es_se)),
#                 width = 0) +
#   geom_errorbar(aes(xmin = psi - 1.96*psi_sd,
#                     xmax = psi + 1.96*psi_sd),
#                 width = 0) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         aspect.ratio = 1) +
#   coord_cartesian(ylim = c(min(exp(antisense_oligos$del_es[-which(antisense_oligos$outside_exon)])),
#                            max(exp(antisense_oligos$del_es[-which(antisense_oligos$outside_exon)]))),
#                   xlim = c(min(antisense_oligos$psi[-which(antisense_oligos$outside_exon)]),
#                            max(antisense_oligos$psi[-which(antisense_oligos$outside_exon)]))
#   ) +
#   scale_colour_manual(values = c("gray20", "firebrick2")) +
#   ylab("deletion enrichment score (exp)") +
#   xlab("antisense oligonucleotide PSI") +
#   ggtitle("rho = 0.7531447; p = 0.01914")
# 
# 
# ggsave("126_AON_enrichment_score.pdf",
#        height = 6,
#        width = 6,
#        useDingbats = F)




ggplot(data = antisense_oligos %>%
         filter(! outside_exon),
       mapping = aes(x = psi,
                     y = del_psi)) +
  geom_smooth(data = subset(antisense_oligos,
                            outside_exon == F),
              method = "lm",
              colour = "black") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = min_del_psi,
                    ymax = max_del_psi),
                width = 0) +
  geom_errorbar(aes(xmin = psi - 1.96*psi_sd,
                    xmax = psi + 1.96*psi_sd),
                width = 0) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  coord_cartesian(ylim = c(min(antisense_oligos$del_psi[-which(antisense_oligos$outside_exon)]),
                           max(antisense_oligos$del_psi[-which(antisense_oligos$outside_exon)])),
                  xlim = c(min(antisense_oligos$psi[-which(antisense_oligos$outside_exon)]),
                           max(antisense_oligos$psi[-which(antisense_oligos$outside_exon)]))
  ) +
  scale_colour_manual(values = c("gray20", "firebrick2")) +
  ylab("deletion PSI") +
  xlab("antisense oligonucleotide PSI") +
  ggtitle("rho = 0.7531447; p = 0.01914")


ggsave("126B_AON_PSI.pdf",
       height = 6,
       width = 6,
       useDingbats = F)





ggplot(data = antisense_oligos %>%
         filter(! outside_exon),
       mapping = aes(x = psi,
                     y = del_spliceAI)) +
  geom_smooth(data = subset(antisense_oligos,
                            outside_exon == F),
              method = "lm",
              colour = "black") +
  geom_point(size = 2) +
  geom_errorbar(aes(xmin = psi - 1.96*psi_sd,
                    xmax = psi + 1.96*psi_sd),
                width = 0) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  coord_cartesian(ylim = c(min(antisense_oligos$del_spliceAI[-which(antisense_oligos$outside_exon)]),
                           max(antisense_oligos$del_spliceAI[-which(antisense_oligos$outside_exon)])),
                  xlim = c(min(antisense_oligos$psi[-which(antisense_oligos$outside_exon)]),
                           max(antisense_oligos$psi[-which(antisense_oligos$outside_exon)]))
  ) +
  scale_colour_manual(values = c("gray20", "firebrick2")) +
  ylab("deletion effect (SpliceAI score)") +
  xlab("antisense oligonucleotide PSI") +
  ggtitle("rho = 0.7352941; p = 0.02398")


ggsave("126C_AON_SPLICEAI.pdf",
       height = 6,
       width = 6,
       useDingbats = F)
