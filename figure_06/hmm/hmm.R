
setwd("/omics/groups/OE0433/internal/pablo/projects/sandbox/illumina/")

library(tidyverse)
library(depmixS4)

dpsi2A <- function(dpsi, wtpsi) {

  dpsi <- dpsi*100
  
  if (dpsi > 0) {
    if ((100 - dpsi) < wtpsi) {
      dpsi <- 99.99 - wtpsi
    }
  }
  
  else if (dpsi < 0) {
    if (wtpsi + dpsi < 0) {
      dpsi <- 0.01 - wtpsi
    }
  }
  
  ((wtpsi^2) - (100*wtpsi) + (dpsi*wtpsi) - (100*dpsi)) / (wtpsi*(dpsi + wtpsi - 100))
}


load("001_master_bed_table_exons_in_80pcnt_samples_with_refseqs_and_psiVals.RData")




# only look at exons of length 50-150
# which(grepl(pattern = "FAS", x = bed_file_exons_in_80_pcnt_of_samples$exon_id) & bed_file_exons_in_80_pcnt_of_samples$length == 63)
# 
# idx <- c(11745, which(bed_file_exons_in_80_pcnt_of_samples$length == 100))
idx <- which(bed_file_exons_in_80_pcnt_of_samples$length >= 50 & bed_file_exons_in_80_pcnt_of_samples$length <= 200)
bed_file_exons_in_80_pcnt_of_samples <- bed_file_exons_in_80_pcnt_of_samples[idx,]

my_exons <- bed_file_exons_in_80_pcnt_of_samples$exon_id









list_of_timeseries_spliceAI <- vector(mode = "list", length = length(my_exons))
list_of_timeseries_logA <- vector(mode = "list", length = length(my_exons))
wt_psi_vector <- c()






for (i in 1:length(my_exons)) {
  
  # i <- 7423
  print(i)
  
  
  # what exon is this?
  this_exon <- my_exons[i]

  # check that file exists
  spliceAI_predictions_file <- paste("data/exon_predictions/",
                                     gsub(pattern = "/", replacement = "-", x = this_exon),
                                     ".tsv",
                                     sep = "",
                                     collapse = "")
  
  # only continue if this file exists (i.e. for which spliceAI provided predictions)
  if (! file.exists(spliceAI_predictions_file)) {
    wt_psi_vector <- c(wt_psi_vector,
                       NA)
    list_of_timeseries_spliceAI[[i]] <- NA
    
    list_of_timeseries_logA[[i]] <- NA
    next
  }
  
  # will now get the mean psi for this exon
  bed_idx <- which(bed_file_exons_in_80_pcnt_of_samples$exon_id == this_exon)
  this_exon_psi <- 100 * mean(x = suppressWarnings(as.numeric(bed_file_exons_in_80_pcnt_of_samples[bed_idx,9:ncol(bed_file_exons_in_80_pcnt_of_samples)])), na.rm = TRUE)
  wt_psi_vector <- c(wt_psi_vector,
                     this_exon_psi)
  
  
  # load spliceAI predictions
  spliceAI_predictions_tibble <- read_delim(file = spliceAI_predictions_file,
                                            delim = "\t",
                                            col_types = cols(
                                              exon = col_character(),
                                              exon_length = col_double(),
                                              mutation = col_character(),
                                              prediction = col_double()
                                            ))
  
  # subset to only k4 deletions (the largest deletions provided by spliceAI with the pre-computed scores)
  spliceAI_predictions_tibble <- spliceAI_predictions_tibble %>%
    filter(grepl(pattern = "Del_k4",
                 x = mutation))
  
  
  mutation_positions <- sapply(X = spliceAI_predictions_tibble$mutation,
                               FUN = function(x){
                                 start_position <- as.numeric(strsplit(x = strsplit(x = x, split = "_")[[1]][3], split = "-")[[1]][1])
                                 k <- 4
                                 start_position + (k/2) - 0.5
                               })
  
  spliceAI_predictions <- spliceAI_predictions_tibble$prediction[order(mutation_positions)]
  
  list_of_timeseries_spliceAI[[i]] <- spliceAI_predictions
  list_of_timeseries_logA[[i]] <- log2(sapply(X = spliceAI_predictions,
                                              FUN = dpsi2A,
                                              wtpsi = this_exon_psi)) - log2(dpsi2A(dpsi = 0, wtpsi = this_exon_psi))

}



names(list_of_timeseries_logA) <- my_exons
list_of_timeseries_logA <- list_of_timeseries_logA[-which(is.na(list_of_timeseries_logA))]

names(list_of_timeseries_spliceAI) <- my_exons
list_of_timeseries_spliceAI <- list_of_timeseries_spliceAI[-which(is.na(list_of_timeseries_spliceAI))]






names(wt_psi_vector) <- my_exons
wt_psi_vector <- wt_psi_vector[-which(is.na(wt_psi_vector))]


my_timeseries_vector <- unlist(list_of_timeseries_spliceAI)
# my_timeseries_vector <- unlist(list_of_timeseries_logA)

timeseries_lengths <- sapply(X = list_of_timeseries_logA,
                             FUN = length)
hmm_model <- depmix(my_timeseries_vector ~ 1,
                    nstates = 3,
                    ntimes = timeseries_lengths)

set.seed(123)
hmm_fit <- fit(hmm_model)

summary(hmm_fit)
getpars(hmm_fit)

model_parameters <- unlist(getpars(hmm_fit))
model_parameters[17] <- 0.15
model_parameters[15] <- 0
model_parameters[13] <- -0.15

#sd's
model_parameters[18] <- 0.1
model_parameters[16] <- 0.025
model_parameters[14] <- 0.1

conpat <- rep(1, 18)
conpat[17] <- 0
conpat[15] <- 0
conpat[13] <- 0

conpat[18] <- 0
conpat[16] <- 0
conpat[14] <- 0

hmm_fit <- setpars(hmm_model, model_parameters)

hmm_fit2 <- fit(hmm_fit, equal = conpat)

save(hmm_fit2, file = "hmm_fit2.RData")

load("hmm_fit2.RData")
# summary(hmm_fit2)



# which(grepl(pattern = "FAS", x = bed_file_exons_in_80_pcnt_of_samples$exon_id) & bed_file_exons_in_80_pcnt_of_samples$length == 63)
# bed_file_exons_in_80_pcnt_of_samples$exon_id[7423]
# plot(list_of_timeseries_logA[["SE_10_+_90770357_90770510_90770572_90771756_FAS"]], type = "l")





hidden_states_split_by_exon <- split(x = hmm_fit2@posterior$state,
                                     f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))
names(hidden_states_split_by_exon) <- names(list_of_timeseries_spliceAI)

save(hidden_states_split_by_exon,
     file = "hmm_hidden_states_split_by_exon.RData")

# hidden_states_split_by_exon[["SE_10_+_90770357_90770510_90770572_90771756_FAS"]]
# 
# plot(list_of_timeseries_spliceAI[["SE_10_+_90770357_90770510_90770572_90771756_FAS"]], type = "l")
# 
# fas_hidden_states <- split(x = hmm_fit2@posterior$state,
#                            f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[[1]]
# plot(x = (1:60)+1.5,
#      y = hidden_states_split_by_exon[["SE_10_+_90770357_90770510_90770572_90771756_FAS"]],
#      type = "l",
#      ylab = "hidden state",
#      xlab = "exon position")
# 
# 
# 
# 
# plot(list_of_timeseries_spliceAI[["SE_10_+_90770357_90770510_90770572_90771756_FAS"]])

plot_tibble <- tibble(position = (1:60)+1.5,
                      state = hidden_states_split_by_exon[["SE_10_+_90770357_90770510_90770572_90771756_FAS"]])

ggplot(data = plot_tibble,
       mapping = aes(x = position,
                     y = state)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 20/63) +
  scale_x_continuous(limits = c(2,62),
                     expand = c(0.01, 0.01),
                     breaks = 1:63)


ggsave(filename = "hmm_fas_exon_modules_plot.pdf",
       width = 10,
       height = 10,
       useDingbats = F)








# 
# 
# mean(hidden_states_split_by_exon[["SE_10_+_90770357_90770510_90770572_90771756_FAS"]] == 2)
ggplot(data = tibble(psi = wt_psi_vector,
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths)),
                                      FUN = function(x){
                                        100 * mean(x == 2)
                                      })),
       mapping = aes(y = percent,
                     x = psi < 90,
                     fill = psi < 90)) +
  geom_violin(alpha = 0.5) +
  geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("gray80", "#F9C45F"))


# 
# ggsave(filename = "hmm_percent_sequence_N.pdf",
#        width = 6,
#        height = 6,
#        useDingbats = F)


# 
# plot(y = sapply(X = split(x = hmm_fit2@posterior$state,
#                           f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths)),
#                 FUN = function(x){
#                   100 * mean(x == 2)
#                 }),
#      x = wt_psi_vector)


exon_lengths <- sapply(X = split(x = hmm_fit2@posterior$state,
                                 f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths)),
                       FUN = function(x){
                         
                         exon_length <- length(x) + 3
                         
                         exon_length
                       })
save(exon_lengths,
     file = "hmm_exon_lengths.RData")

fraction_n_state_tibble <- tibble(psi = wt_psi_vector,
                                  percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                             f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths)),
                                                   FUN = function(x){
                                                     100 * mean(x == 2)
                                                   }),
                                  length = exon_lengths)

median(fraction_n_state_tibble$percent[fraction_n_state_tibble$psi<90 & fraction_n_state_tibble$length <= 100])
median(fraction_n_state_tibble$percent[fraction_n_state_tibble$psi>=90 & fraction_n_state_tibble$length <= 100])


ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths <= 100],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths <= 100],
                                      FUN = function(x){
                                        100 * mean(x == 2)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +
  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_N_violin_plots_under100.pdf",
       width = 8,
       height = 8,
       useDingbats = F)








ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths <= 100],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths <= 100],
                                      FUN = function(x){
                                        100 * mean(x == 3)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +
  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_S_violin_plots_under100.pdf",
       width = 8,
       height = 8,
       useDingbats = F)



ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths <= 100],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths <= 100],
                                      FUN = function(x){
                                        100 * mean(x == 1)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +
  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_E_violin_plots_under100.pdf",
       width = 8,
       height = 8,
       useDingbats = F)











ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 100 & exon_lengths <= 150],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths > 100 & exon_lengths <= 150],
                                      FUN = function(x){
                                        100 * mean(x == 2)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_N_violin_plots_100to150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)





ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 100 & exon_lengths <= 150],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths > 100 & exon_lengths <= 150],
                                      FUN = function(x){
                                        100 * mean(x == 3)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_S_violin_plots_100to150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)




ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 100 & exon_lengths <= 150],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths > 100 & exon_lengths <= 150],
                                      FUN = function(x){
                                        100 * mean(x == 1)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_E_violin_plots_100to150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)







ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 150],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths > 150],
                                      FUN = function(x){
                                        100 * mean(x == 2)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_N_violin_plots_over150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)






ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 150],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths > 150],
                                      FUN = function(x){
                                        100 * mean(x == 3)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_S_violin_plots_over150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)




ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 150],
                     percent = sapply(X = split(x = hmm_fit2@posterior$state,
                                                f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths))[exon_lengths > 150],
                                      FUN = function(x){
                                        100 * mean(x == 1)
                                      })),
       mapping = aes(y = percent,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_percent_sequence_E_violin_plots_over150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)




fas_exon_length <- 63
fas_sre_lengths <- rle(fas_hidden_states)$lengths
fas_changes_of_state <- length(rle(fas_hidden_states)$values)-1


fas_changes_of_state/fas_exon_length



# 
# # changes of state per nucleotide
# changes_of_state <- sapply(X = split(x = hmm_fit2@posterior$state,
#                                      f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths)),
#                            FUN = function(x){
#                              
#                              exon_length <- length(x) + 3
#                              
#                              rle_object <- rle(x = x)
#                              changes_of_state <- length(rle_object$values) - 1
#                              
#                              
#                              changes_of_state/exon_length
#                            })



# changes of state per nucleotide
changes_of_state <- sapply(X = split(x = hmm_fit2@posterior$state,
                                     f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths)),
                           FUN = function(x){
                             
                             exon_length <- length(x) + 3
                             
                             rle_object <- rle(x = x)
                             
                             
                             temp_tibble <- tibble(value = rle_object$values,
                                                   next_value = c(rle_object$values[-1],NA))
                             
                             
                             one_to_three <- as.numeric(table(temp_tibble$next_value[temp_tibble$value == 1])["3"])
                             three_to_one <- as.numeric(table(temp_tibble$next_value[temp_tibble$value == 3])["1"])
                             
                             if (is.na(one_to_three)) {
                               one_to_three <- 0
                             }
                             
                             if (is.na(three_to_one)) {
                               three_to_one <- 0
                             }
                             
                             
                             changes_of_state <- one_to_three + three_to_one
                             
                             
                             changes_of_state/exon_length
                           })



# 
# ggplot(data = tibble(psi = wt_psi_vector,
#                      changes = changes_of_state),
#        mapping = aes(x = changes,
#                      fill = psi < 90)) +
#   geom_density(alpha = 0.5) +
#   geom_vline(xintercept = changes_of_state[1]) +
#   theme_bw() +
#   theme(aspect.ratio = 1,
#         panel.grid = element_blank()) +
#   scale_fill_manual(values = c("gray90", "#F9C45F")) +
#   xlab("changes of state per nt")
# 
# 
# ggsave(filename = "hmm_changes_state_per_nt.pdf",
#        width = 6,
#        height = 6,
#        useDingbats = F)




ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths <= 100],
                     changes_per_nt = changes_of_state[exon_lengths <= 100]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_ES_states_per_nt_violin_plots_under100.pdf",
       width = 8,
       height = 8,
       useDingbats = F)



ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 100 & exon_lengths <= 150],
                     changes_per_nt = changes_of_state[exon_lengths > 100 & exon_lengths <= 150]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_ES_states_per_nt_violin_plots_100_to_150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)


ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 150],
                     changes_per_nt = changes_of_state[exon_lengths > 150]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_ES_states_per_nt_violin_plots_over_150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)


changes_state_tibble <- tibble(psi = wt_psi_vector,
                               changes_per_nt = changes_of_state,
                               length = exon_lengths)


median(changes_state_tibble$changes_per_nt[changes_state_tibble$psi<90 & changes_state_tibble$length <= 100])
median(changes_state_tibble$changes_per_nt[changes_state_tibble$psi>=90 & changes_state_tibble$length <= 100])

hist(changes_state_tibble$changes_per_nt[changes_state_tibble$psi<90 & changes_state_tibble$length <= 100])
hist(changes_state_tibble$changes_per_nt[changes_state_tibble$psi>=90 & changes_state_tibble$length <= 100])


wilcox.test(x = changes_state_tibble$changes_per_nt[changes_state_tibble$psi<90 & changes_state_tibble$length <= 100],
            y = changes_state_tibble$changes_per_nt[changes_state_tibble$psi>=90 & changes_state_tibble$length <= 100])













# changes of state per nucleotide
changes_of_EN_state <- sapply(X = split(x = hmm_fit2@posterior$state,
                                     f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths)),
                           FUN = function(x){
                             
                             exon_length <- length(x) + 3
                             
                             rle_object <- rle(x = x)
                             
                             
                             temp_tibble <- tibble(value = rle_object$values,
                                                   next_value = c(rle_object$values[-1],NA))
                             
                             
                             two_to_one <- as.numeric(table(temp_tibble$next_value[temp_tibble$value == 2])["1"])
                             one_to_two <- as.numeric(table(temp_tibble$next_value[temp_tibble$value == 1])["2"])
                             
                             if (is.na(two_to_one)) {
                               two_to_one <- 0
                             }
                             
                             if (is.na(one_to_two)) {
                               one_to_two <- 0
                             }
                             
                             
                             changes_of_state <- two_to_one + one_to_two
                             
                             
                             changes_of_state/exon_length
                           })






ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths <= 100],
                     changes_per_nt = changes_of_EN_state[exon_lengths <= 100]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_EN_states_per_nt_violin_plots_under100.pdf",
       width = 8,
       height = 8,
       useDingbats = F)



ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 100 & exon_lengths <= 150],
                     changes_per_nt = changes_of_EN_state[exon_lengths > 100 & exon_lengths <= 150]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_EN_states_per_nt_violin_plots_100_to_150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)


ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 150],
                     changes_per_nt = changes_of_EN_state[exon_lengths > 150]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_EN_states_per_nt_violin_plots_over_150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)





changes_EN_state_tibble <- tibble(psi = wt_psi_vector,
                                  changes_per_nt = changes_of_EN_state,
                                  length = exon_lengths)

median(changes_EN_state_tibble$changes_per_nt[changes_EN_state_tibble$psi<90 & changes_EN_state_tibble$length <= 100])
median(changes_EN_state_tibble$changes_per_nt[changes_EN_state_tibble$psi>=90 & changes_EN_state_tibble$length <= 100])

hist(changes_state_tibble$changes_per_nt[changes_state_tibble$psi<90 & changes_state_tibble$length <= 100])
hist(changes_state_tibble$changes_per_nt[changes_state_tibble$psi>=90 & changes_state_tibble$length <= 100])


wilcox.test(x = changes_EN_state_tibble$changes_per_nt[changes_EN_state_tibble$psi<90 & changes_EN_state_tibble$length <= 100],
            y = changes_EN_state_tibble$changes_per_nt[changes_EN_state_tibble$psi>=90 & changes_EN_state_tibble$length <= 100])



median(changes_EN_state_tibble$changes_per_nt[changes_EN_state_tibble$psi<90 & changes_EN_state_tibble$length > 100])
median(changes_EN_state_tibble$changes_per_nt[changes_EN_state_tibble$psi>=90 & changes_EN_state_tibble$length > 100])

wilcox.test(x = changes_EN_state_tibble$changes_per_nt[changes_EN_state_tibble$psi<90 & changes_EN_state_tibble$length > 100],
            y = changes_EN_state_tibble$changes_per_nt[changes_EN_state_tibble$psi>=90 & changes_EN_state_tibble$length > 100])















# changes of state per nucleotide
changes_of_SN_state <- sapply(X = split(x = hmm_fit2@posterior$state,
                                        f = rep(1:length(list_of_timeseries_spliceAI), timeseries_lengths)),
                              FUN = function(x){
                                
                                exon_length <- length(x) + 3
                                
                                rle_object <- rle(x = x)
                                
                                
                                temp_tibble <- tibble(value = rle_object$values,
                                                      next_value = c(rle_object$values[-1],NA))
                                
                                
                                two_to_thr <- as.numeric(table(temp_tibble$next_value[temp_tibble$value == 2])["3"])
                                thr_to_two <- as.numeric(table(temp_tibble$next_value[temp_tibble$value == 3])["2"])
                                
                                if (is.na(two_to_thr)) {
                                  two_to_thr <- 0
                                }
                                
                                if (is.na(thr_to_two)) {
                                  thr_to_two <- 0
                                }
                                
                                
                                changes_of_state <- two_to_thr + thr_to_two
                                
                                
                                changes_of_state/exon_length
                              })



ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths <= 100],
                     changes_per_nt = changes_of_SN_state[exon_lengths <= 100]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_SN_states_per_nt_violin_plots_under100.pdf",
       width = 8,
       height = 8,
       useDingbats = F)



ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 100 & exon_lengths <= 150],
                     changes_per_nt = changes_of_SN_state[exon_lengths > 100 & exon_lengths <= 150]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_SN_states_per_nt_violin_plots_100_to_150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)


ggplot(data = tibble(psi = as.factor(findInterval(x = wt_psi_vector,
                                                  vec = seq(0,100,10), all.inside = T))[exon_lengths > 150],
                     changes_per_nt = changes_of_SN_state[exon_lengths > 150]),
       mapping = aes(y = changes_per_nt,
                     x = psi,
                     group = psi)) +
  geom_violin(alpha = 0.5, scale = "width", fill = "gray80") +
  # geom_boxplot(notch = T, width = 0.2, outlier.shape = NA, coef = 0) +
  stat_summary(geom = "point",
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +  # geom_vline(xintercept = 16.66667) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(breaks=c("1","2","3","4","5", "6", "7", "8", "9", "10"),
                   labels=c("[0 - 10)", "[10 - 20)", "[20 - 30)", "[30 - 40)", "[40 - 50]",
                            "[50 - 60)", "[60 - 70)", "[70 - 80)", "[80 - 90)", "[90 - 100]"))

ggsave(filename = "hmm_changes_SN_states_per_nt_violin_plots_over_150.pdf",
       width = 8,
       height = 8,
       useDingbats = F)




changes_SN_state_tibble <- tibble(psi = wt_psi_vector,
                                  changes_per_nt = changes_of_SN_state,
                                  length = exon_lengths)

median(changes_SN_state_tibble$changes_per_nt[changes_SN_state_tibble$psi<90 & changes_SN_state_tibble$length <= 100])
median(changes_SN_state_tibble$changes_per_nt[changes_SN_state_tibble$psi>=90 & changes_SN_state_tibble$length <= 100])

hist(changes_state_tibble$changes_per_nt[changes_state_tibble$psi<90 & changes_state_tibble$length <= 100])
hist(changes_state_tibble$changes_per_nt[changes_state_tibble$psi>=90 & changes_state_tibble$length <= 100])


wilcox.test(x = changes_SN_state_tibble$changes_per_nt[changes_SN_state_tibble$psi<90 & changes_SN_state_tibble$length <= 100],
            y = changes_SN_state_tibble$changes_per_nt[changes_SN_state_tibble$psi>=90 & changes_SN_state_tibble$length <= 100])



median(changes_SN_state_tibble$changes_per_nt[changes_SN_state_tibble$psi<90 & changes_SN_state_tibble$length > 100])
median(changes_SN_state_tibble$changes_per_nt[changes_SN_state_tibble$psi>=90 & changes_SN_state_tibble$length > 100])

wilcox.test(x = changes_SN_state_tibble$changes_per_nt[changes_SN_state_tibble$psi<90 & changes_SN_state_tibble$length > 100],
            y = changes_SN_state_tibble$changes_per_nt[changes_SN_state_tibble$psi>=90 & changes_SN_state_tibble$length > 100])


























sre_lengths <- unlist(sapply(X = split(x = hmm_fit2@posterior$state,
                                       f = rep(1:length(list_of_timeseries_logA), timeseries_lengths)),
                             FUN = function(x){
                               rle_object <- rle(x = x)
                               
                               idx_enhancers <- which(rle_object$values == 1)
                               idx_silencers <- which(rle_object$values == 3)
                               
                               enhancer_lengths <- rle_object$lengths[idx_enhancers]
                               silencer_lengths <- rle_object$lengths[idx_silencers]
                               
                               c(enhancer_lengths, silencer_lengths)
                             }))

ggplot(data = tibble(length = sre_lengths),
       mapping = aes(x = length)) +
  geom_density(fill = "gray90", alpha = 0.5, bw = 2) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  xlab("length (nt)")

ggsave(filename = "hmm_sre_lengths.pdf",
       width = 8,
       height = 8,
       useDingbats = F)
median(sre_lengths)
mean(sre_lengths)



enh_lengths <- unlist(sapply(X = split(x = hmm_fit2@posterior$state,
                                       f = rep(1:length(list_of_timeseries_logA), timeseries_lengths)),
                             FUN = function(x){
                               rle_object <- rle(x = x)
                               
                               idx_enhancers <- which(rle_object$values == 1)
                               # idx_silencers <- which(rle_object$values == 3)
                               
                               enhancer_lengths <- rle_object$lengths[idx_enhancers]
                               # silencer_lengths <- rle_object$lengths[idx_silencers]
                               
                               enhancer_lengths
                             }))

ggplot(data = tibble(length = enh_lengths),
       mapping = aes(x = length)) +
  geom_density(fill = "gray90", alpha = 0.5, bw = 2) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  xlab("length (nt)")

ggsave(filename = "hmm_enh_lengths.pdf",
       width = 8,
       height = 8,
       useDingbats = F)
median(enh_lengths)
mean(enh_lengths)




sil_lengths <- unlist(sapply(X = split(x = hmm_fit2@posterior$state,
                                       f = rep(1:length(list_of_timeseries_logA), timeseries_lengths)),
                             FUN = function(x){
                               rle_object <- rle(x = x)
                               
                               # idx_enhancers <- which(rle_object$values == 1)
                               idx_silencers <- which(rle_object$values == 3)
                               
                               # enhancer_lengths <- rle_object$lengths[idx_enhancers]
                               silencer_lengths <- rle_object$lengths[idx_silencers]
                               
                               silencer_lengths
                             }))

ggplot(data = tibble(length = sil_lengths),
       mapping = aes(x = length)) +
  geom_density(fill = "gray90", alpha = 0.5, bw = 2) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  xlab("length (nt)")

ggsave(filename = "hmm_sil_lengths.pdf",
       width = 8,
       height = 8,
       useDingbats = F)
median(sil_lengths)
mean(sil_lengths)



ggplot(data = tibble(length = c(sil_lengths, enh_lengths),
                     class = c(rep("S", length(sil_lengths)),
                               rep("E", length(enh_lengths)))),
       mapping = aes(x = length,
                     colour = class)) +
  geom_density(bw = 1) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank()) +
  xlab("length (nt)")


ggsave(filename = "hmm_enh_sil_lengths.pdf",
       width = 8,
       height = 8,
       useDingbats = F)


rep(c("a", "b", "c"), c(1,2,3))


unique(plot_tibble$psi_group)

plot_tibble <- tibble(spliceAI = abs(unlist(list_of_timeseries_spliceAI)),
                      psi_group = rep(findInterval(x = wt_psi_vector,
                                                   vec = seq(0,100,10), all.inside = TRUE),
                                      sapply(list_of_timeseries_spliceAI, length)))

ggplot(data = plot_tibble,
       mapping = aes(x = spliceAI,
                     colour = as.factor(psi_group))) +
  geom_density(bw = 0.01) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank()) +
  scale_colour_manual(values = c("blue", "gray55", "gray60", "gray65", "gray70",
                                 "gray75", "gray80", "gray85", "gray90", "red"))

ggsave(filename = "hmm_spliceAI_effects.pdf",
       width = 8,
       height = 8,
       useDingbats = F)

plot(density(unlist(list_of_timeseries_spliceAI)))




























###### State per position


extract_elements <- function(list_of_vectors, position = 1) {
  extracted_elements <- sapply(list_of_vectors, function(vec) {
    if (position > 0) {
      vec[position]
    } else {
      vec[length(vec) + position + 1]
    }
  })
  return(extracted_elements)
}


elements_first_20 <- sapply(X = 1:20,
                            FUN = function(x){
                              extract_elements(list_of_vectors = hidden_states_split_by_exon,
                                               position = x)
                            })

elements_last_20 <- sapply(X = -20:-1,
                           FUN = function(x){
                             extract_elements(list_of_vectors = hidden_states_split_by_exon,
                                              position = x)
                           })




percentages_first_20 <- apply(X = apply(X = elements_first_20,
                                        MARGIN = 2,
                                        FUN = table),
                              MARGIN = 2,
                              FUN = function(y) {
                                100 * y / sum(y)
                              })

percentages_last_20 <- apply(X = apply(X = elements_last_20,
                                       MARGIN = 2,
                                       FUN = table),
                             MARGIN = 2,
                             FUN = function(y) {
                               100 * y / sum(y)
                             })





tibble_freq_first_20 <- tibble(percent = c(percentages_first_20),
                               position = rep(1:20, each = 3),
                               state = rep(1:3, 20))



tibble_freq_last_20 <- tibble(percent = c(percentages_last_20),
                              position = rep(-20:-1, each = 3),
                              state = rep(1:3, 20))


ggplot(data = tibble_freq_first_20,
       mapping = aes(x = factor(position),
                     y = percent,
                     fill = factor(state,
                                   levels = c("3", "2", "1")))) +
  geom_bar(stat = "identity") +
  labs(title = "States in the first 20 4mers of all exons",
       x = "4mer window",
       y = "percent") +
  scale_fill_discrete(name = "state",
                      labels = c("S", "N", "E")) +
  scale_x_discrete(labels = sapply(X = 1:20,
                                   FUN = function(x){
                                     paste(x, " - ", x+3, sep = "", collapse = "")
                                   }), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 0.65)

ggsave(filename = "hmm_first_20_positions.pdf",
       width = 6,
       height = 6,
       useDingbats = F)


ggplot(data = tibble_freq_last_20,
       mapping = aes(x = factor(position),
                     y = percent,
                     fill = factor(state,
                                   levels = c("3", "2", "1")))) +
  geom_bar(stat = "identity") +
  labs(title = "States in the last 20 4mers of all exons",
       x = "4mer window",
       y = "percent") +
  scale_fill_discrete(name = "state",
                      labels = c("S", "N", "E")) +
  scale_x_discrete(labels = sapply(X = 1:20,
                                   FUN = function(x){
                                     paste(x, " - ", x+3, sep = "", collapse = "")
                                   }), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 0.65)

ggsave(filename = "hmm_last_20_positions.pdf",
       width = 6,
       height = 6,
       useDingbats = F)




























# SHORT EXONS



extract_elements <- function(list_of_vectors, position = 1) {
  extracted_elements <- sapply(list_of_vectors, function(vec) {
    if (position > 0) {
      vec[position]
    } else {
      vec[length(vec) + position + 1]
    }
  })
  return(extracted_elements)
}


elements_first_20 <- sapply(X = 1:20,
                            FUN = function(x){
                              extract_elements(list_of_vectors = hidden_states_split_by_exon[which(exon_lengths < 100)],
                                               position = x)
                            })

elements_last_20 <- sapply(X = -20:-1,
                           FUN = function(x){
                             extract_elements(list_of_vectors = hidden_states_split_by_exon[which(exon_lengths < 100)],
                                              position = x)
                           })




percentages_first_20 <- apply(X = apply(X = elements_first_20,
                                        MARGIN = 2,
                                        FUN = table),
                              MARGIN = 2,
                              FUN = function(y) {
                                100 * y / sum(y)
                              })

percentages_last_20 <- apply(X = apply(X = elements_last_20,
                                       MARGIN = 2,
                                       FUN = table),
                             MARGIN = 2,
                             FUN = function(y) {
                               100 * y / sum(y)
                             })





tibble_freq_first_20 <- tibble(percent = c(percentages_first_20),
                               position = rep(1:20, each = 3),
                               state = rep(1:3, 20))



tibble_freq_last_20 <- tibble(percent = c(percentages_last_20),
                              position = rep(-20:-1, each = 3),
                              state = rep(1:3, 20))


ggplot(data = tibble_freq_first_20,
       mapping = aes(x = factor(position),
                     y = percent,
                     fill = factor(state,
                                   levels = c("3", "2", "1")))) +
  geom_bar(stat = "identity") +
  labs(title = "States in the first 20 4mers of all exons",
       x = "4mer window",
       y = "percent") +
  scale_fill_discrete(name = "state",
                      labels = c("S", "N", "E")) +
  scale_x_discrete(labels = sapply(X = 1:20,
                                   FUN = function(x){
                                     paste(x, " - ", x+3, sep = "", collapse = "")
                                   }), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 0.65)

ggsave(filename = "hmm_first_20_positions_UNDER100.pdf",
       width = 6,
       height = 6,
       useDingbats = F)


ggplot(data = tibble_freq_last_20,
       mapping = aes(x = factor(position),
                     y = percent,
                     fill = factor(state,
                                   levels = c("3", "2", "1")))) +
  geom_bar(stat = "identity") +
  labs(title = "States in the last 20 4mers of all exons",
       x = "4mer window",
       y = "percent") +
  scale_fill_discrete(name = "state",
                      labels = c("S", "N", "E")) +
  scale_x_discrete(labels = sapply(X = 1:20,
                                   FUN = function(x){
                                     paste(x, " - ", x+3, sep = "", collapse = "")
                                   }), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 0.65)

ggsave(filename = "hmm_last_20_positions_UNDER100.pdf",
       width = 6,
       height = 6,
       useDingbats = F)




























# MEDIUM EXONS



extract_elements <- function(list_of_vectors, position = 1) {
  extracted_elements <- sapply(list_of_vectors, function(vec) {
    if (position > 0) {
      vec[position]
    } else {
      vec[length(vec) + position + 1]
    }
  })
  return(extracted_elements)
}


elements_first_20 <- sapply(X = 1:20,
                            FUN = function(x){
                              extract_elements(list_of_vectors = hidden_states_split_by_exon[which(exon_lengths < 150 & exon_lengths >= 100)],
                                               position = x)
                            })

elements_last_20 <- sapply(X = -20:-1,
                           FUN = function(x){
                             extract_elements(list_of_vectors = hidden_states_split_by_exon[which(exon_lengths < 150 & exon_lengths >= 100)],
                                              position = x)
                           })




percentages_first_20 <- apply(X = apply(X = elements_first_20,
                                        MARGIN = 2,
                                        FUN = table),
                              MARGIN = 2,
                              FUN = function(y) {
                                100 * y / sum(y)
                              })

percentages_last_20 <- apply(X = apply(X = elements_last_20,
                                       MARGIN = 2,
                                       FUN = table),
                             MARGIN = 2,
                             FUN = function(y) {
                               100 * y / sum(y)
                             })





tibble_freq_first_20 <- tibble(percent = c(percentages_first_20),
                               position = rep(1:20, each = 3),
                               state = rep(1:3, 20))



tibble_freq_last_20 <- tibble(percent = c(percentages_last_20),
                              position = rep(-20:-1, each = 3),
                              state = rep(1:3, 20))


ggplot(data = tibble_freq_first_20,
       mapping = aes(x = factor(position),
                     y = percent,
                     fill = factor(state,
                                   levels = c("3", "2", "1")))) +
  geom_bar(stat = "identity") +
  labs(title = "States in the first 20 4mers of all exons",
       x = "4mer window",
       y = "percent") +
  scale_fill_discrete(name = "state",
                      labels = c("S", "N", "E")) +
  scale_x_discrete(labels = sapply(X = 1:20,
                                   FUN = function(x){
                                     paste(x, " - ", x+3, sep = "", collapse = "")
                                   }), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 0.65)

ggsave(filename = "hmm_first_20_positions_100to150.pdf",
       width = 6,
       height = 6,
       useDingbats = F)


ggplot(data = tibble_freq_last_20,
       mapping = aes(x = factor(position),
                     y = percent,
                     fill = factor(state,
                                   levels = c("3", "2", "1")))) +
  geom_bar(stat = "identity") +
  labs(title = "States in the last 20 4mers of all exons",
       x = "4mer window",
       y = "percent") +
  scale_fill_discrete(name = "state",
                      labels = c("S", "N", "E")) +
  scale_x_discrete(labels = sapply(X = 1:20,
                                   FUN = function(x){
                                     paste(x, " - ", x+3, sep = "", collapse = "")
                                   }), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 0.65)

ggsave(filename = "hmm_last_20_positions_100to150.pdf",
       width = 6,
       height = 6,
       useDingbats = F)
























# LONG EXONS



extract_elements <- function(list_of_vectors, position = 1) {
  extracted_elements <- sapply(list_of_vectors, function(vec) {
    if (position > 0) {
      vec[position]
    } else {
      vec[length(vec) + position + 1]
    }
  })
  return(extracted_elements)
}


elements_first_20 <- sapply(X = 1:20,
                            FUN = function(x){
                              extract_elements(list_of_vectors = hidden_states_split_by_exon[which(exon_lengths >= 150)],
                                               position = x)
                            })

elements_last_20 <- sapply(X = -20:-1,
                           FUN = function(x){
                             extract_elements(list_of_vectors = hidden_states_split_by_exon[which(exon_lengths >= 150)],
                                              position = x)
                           })




percentages_first_20 <- apply(X = apply(X = elements_first_20,
                                        MARGIN = 2,
                                        FUN = table),
                              MARGIN = 2,
                              FUN = function(y) {
                                100 * y / sum(y)
                              })

percentages_last_20 <- apply(X = apply(X = elements_last_20,
                                       MARGIN = 2,
                                       FUN = table),
                             MARGIN = 2,
                             FUN = function(y) {
                               100 * y / sum(y)
                             })





tibble_freq_first_20 <- tibble(percent = c(percentages_first_20),
                               position = rep(1:20, each = 3),
                               state = rep(1:3, 20))



tibble_freq_last_20 <- tibble(percent = c(percentages_last_20),
                              position = rep(-20:-1, each = 3),
                              state = rep(1:3, 20))


ggplot(data = tibble_freq_first_20,
       mapping = aes(x = factor(position),
                     y = percent,
                     fill = factor(state,
                                   levels = c("3", "2", "1")))) +
  geom_bar(stat = "identity") +
  labs(title = "States in the first 20 4mers of all exons",
       x = "4mer window",
       y = "percent") +
  scale_fill_discrete(name = "state",
                      labels = c("S", "N", "E")) +
  scale_x_discrete(labels = sapply(X = 1:20,
                                   FUN = function(x){
                                     paste(x, " - ", x+3, sep = "", collapse = "")
                                   }), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 0.65)

ggsave(filename = "hmm_first_20_positions_OVER150.pdf",
       width = 6,
       height = 6,
       useDingbats = F)


ggplot(data = tibble_freq_last_20,
       mapping = aes(x = factor(position),
                     y = percent,
                     fill = factor(state,
                                   levels = c("3", "2", "1")))) +
  geom_bar(stat = "identity") +
  labs(title = "States in the last 20 4mers of all exons",
       x = "4mer window",
       y = "percent") +
  scale_fill_discrete(name = "state",
                      labels = c("S", "N", "E")) +
  scale_x_discrete(labels = sapply(X = 1:20,
                                   FUN = function(x){
                                     paste(x, " - ", x+3, sep = "", collapse = "")
                                   }), position = "bottom") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        aspect.ratio = 0.65)

ggsave(filename = "hmm_last_20_positions_OVER150.pdf",
       width = 6,
       height = 6,
       useDingbats = F)

















