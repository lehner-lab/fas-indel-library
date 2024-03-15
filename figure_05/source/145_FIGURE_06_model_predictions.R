
library(tidyverse)
library(scales)

setwd("~/CRG/Indel/")



load("model_predictions/SMS_scores/sms_predictions.RData")

plot_tibble <- tibble(id = results_table$id,
                      psi = results_table$psi,
                      prediction = results_table$sms,
                      class = sapply(X = results_table$id,
                                     FUN = function(x){
                                       result <- "?"
                                       
                                       if (grepl(pattern = "Del", x = x)) {
                                         result <- "deletion"
                                       } else if (grepl(pattern = "Ins", x = x)) {
                                         result <- "insertion"
                                       } else if (grepl(pattern = "WT", x = x)) {
                                         result <- "wt"
                                       } else {
                                         result <- "substitution"
                                       }
                                       result
                                     }))

ggplot(mapping = aes(x = prediction,
                     y = psi,
                     colour = class)) +
  geom_point(data = plot_tibble %>%
               filter(class == "insertion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "deletion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "substitution")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  ylab("PSI") +
  xlab("prediction") +
  scale_colour_manual(values = c("#F9C45F", "gray90", "gray30"))

ggsave(filename = "145_sms_scores.pdf",
       height = 7,
       width = 7,
       useDingbats = F)







load("model_predictions/hal/hal_predictions.RData")

plot_tibble <- tibble(id = results_table$id,
                      psi = results_table$psi,
                      prediction = results_table$hal,
                      class = sapply(X = results_table$id,
                                     FUN = function(x){
                                       result <- "?"
                                       
                                       if (grepl(pattern = "Del", x = x)) {
                                         result <- "deletion"
                                       } else if (grepl(pattern = "Ins", x = x)) {
                                         result <- "insertion"
                                       } else if (grepl(pattern = "WT", x = x)) {
                                         result <- "wt"
                                       } else {
                                         result <- "substitution"
                                       }
                                       result
                                     }))

ggplot(mapping = aes(x = prediction,
                     y = psi,
                     colour = class)) +
  geom_point(data = plot_tibble %>%
               filter(class == "insertion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "deletion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "substitution")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  ylab("PSI") +
  xlab("prediction") +
  scale_colour_manual(values = c("#F9C45F", "gray90", "gray30"))

cor.test(x = (plot_tibble %>% filter(class == "deletion"))$psi,
         y = (plot_tibble %>% filter(class == "deletion"))$prediction,
         method = "spearman")

ggsave(filename = "145_hal.pdf",
       height = 7,
       width = 7,
       useDingbats = F)

cor.test(x = (plot_tibble %>% filter(class == "deletion"))$psi,
         y = (plot_tibble %>% filter(class == "deletion"))$prediction,
         method = "spearman")



load("model_predictions/mmsplice/mmsplice_predictions.RData")

plot_tibble <- tibble(id = results_table$id,
                      psi = results_table$psi,
                      prediction = results_table$mmsplice,
                      class = sapply(X = results_table$id,
                                     FUN = function(x){
                                       result <- "?"
                                       
                                       if (grepl(pattern = "Del", x = x)) {
                                         result <- "deletion"
                                       } else if (grepl(pattern = "Ins", x = x)) {
                                         result <- "insertion"
                                       } else if (grepl(pattern = "WT", x = x)) {
                                         result <- "wt"
                                       } else {
                                         result <- "substitution"
                                       }
                                       result
                                     }))

ggplot(mapping = aes(x = prediction,
                     y = psi,
                     colour = class)) +
  geom_point(data = plot_tibble %>%
               filter(class == "insertion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "deletion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "substitution")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  ylab("PSI") +
  xlab("prediction") +
  scale_colour_manual(values = c("#F9C45F", "gray90", "gray30"))

cor.test(x = (plot_tibble %>% filter(class == "deletion"))$psi,
         y = (plot_tibble %>% filter(class == "deletion"))$prediction,
         method = "spearman")

ggsave(filename = "145_mmsplice.pdf",
       height = 7,
       width = 7,
       useDingbats = F)








load("model_predictions/spliceAI/spliceAI_predictions.RData")

plot_tibble <- tibble(id = results_table$id,
                      psi = results_table$psi,
                      prediction = results_table$spliceAI,
                      class = sapply(X = results_table$id,
                                     FUN = function(x){
                                       result <- "?"
                                       
                                       if (grepl(pattern = "Del", x = x)) {
                                         result <- "deletion"
                                       } else if (grepl(pattern = "Ins", x = x)) {
                                         result <- "insertion"
                                       } else if (grepl(pattern = "WT", x = x)) {
                                         result <- "wt"
                                       } else {
                                         result <- "substitution"
                                       }
                                       result
                                     }))

ggplot(mapping = aes(x = prediction,
                     y = psi,
                     colour = class)) +
  geom_point(data = plot_tibble %>%
               filter(class == "insertion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "deletion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "substitution")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  ylab("PSI") +
  xlab("prediction") +
  scale_colour_manual(values = c("#F9C45F", "gray90", "gray30"))


ggsave(filename = "145_spliceAI.pdf",
       height = 7,
       width = 7,
       useDingbats = F)









load("model_predictions/pangolin/pangolin_predictions.RData")

plot_tibble <- tibble(id = results_table$id,
                      psi = results_table$psi,
                      prediction = results_table$pangolin,
                      class = sapply(X = results_table$id,
                                     FUN = function(x){
                                       result <- "?"
                                       
                                       if (grepl(pattern = "Del", x = x)) {
                                         result <- "deletion"
                                       } else if (grepl(pattern = "Ins", x = x)) {
                                         result <- "insertion"
                                       } else if (grepl(pattern = "WT", x = x)) {
                                         result <- "wt"
                                       } else {
                                         result <- "substitution"
                                       }
                                       result
                                     }))

ggplot(mapping = aes(x = prediction,
                     y = psi,
                     colour = class)) +
  geom_point(data = plot_tibble %>%
               filter(class == "insertion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "deletion")) +
  geom_point(data = plot_tibble %>%
               filter(class == "substitution")) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  ylab("PSI") +
  xlab("prediction") +
  scale_colour_manual(values = c("#F9C45F", "gray90", "gray30"))


ggsave(filename = "145_pangolin.pdf",
       height = 7,
       width = 7,
       useDingbats = F)

