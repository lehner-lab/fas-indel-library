
library(tidyverse)
library(scales)
load("data/indels_library.RData")

wt_sequence <- indels_library$sequence[which(indels_library$id == "WT")]
idx_singles <- which(nchar(indels_library$sequence)==63 & indels_library$id != "WT")
density_object <- density(x = indels_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]

indels_library$es <- indels_library$es - library_bias

indels_library$psi <- (49.1/1) * exp(indels_library$es)
indels_library$psi[which(indels_library$psi>100)] <-  100

indels_library$psi_max <- (49.1/1) * exp(indels_library$es + 1.96*indels_library$standard.error)
indels_library$psi_max[which(indels_library$psi_max>100)] <-  100
indels_library$psi_min <- (49.1/1) * exp(indels_library$es - 1.96*indels_library$standard.error)



# 
experimental_validations <- readxl::read_excel(path = "data/Validation_Uber_Indel_clones_140422.xlsx",
                                        skip = 25)


exons_with_new_ss <- c("Indel8", "Indel10", "Indel18", "Indel19", "Indel22")



indels_tibble <- tibble(library_id = sapply(X = experimental_validations$sequence,
                                  FUN = function(x){
                                    idx <- which(indels_library$sequence == x)
                                    
                                    if (length(idx) > 0) {
                                      indels_library$id[idx]
                                    } else {
                                      NA
                                    }
                                    
                                  }),
                      experiment_id = experimental_validations$SampleID,
                      experimental_psi = as.numeric(experimental_validations$PSI),
                      experimental_psi_max = sapply(X = experimental_validations$SampleID,
                                                    FUN = function(x){
                                                      idx <- which(experimental_validations$SampleID == x)
                                                      
                                                      psi <- as.numeric(experimental_validations$PSI)[idx] + 1.96*(as.numeric(experimental_validations$SD)[idx])
                                                      
                                                      if (!is.na(psi)){
                                                        if (psi > 100) {
                                                          psi <- 100
                                                        }
                                                      }
                                                      
                                                      psi
                                                      
                                                    }),
                      experimental_psi_min = sapply(X = experimental_validations$SampleID,
                                                    FUN = function(x){
                                                      idx <- which(experimental_validations$SampleID == x)
                                                      
                                                      psi <- as.numeric(experimental_validations$PSI)[idx] - 1.96*(as.numeric(experimental_validations$SD)[idx])
                                                      
                                                      if(!is.na(psi)){
                                                        if (psi < 0) {
                                                          psi <- 0
                                                        }
                                                      }
                                                      
                                                      psi
                                                      
                                                    }),
                      library_psi = sapply(X = experimental_validations$sequence,
                                           FUN = function(x){
                                             idx <- which(indels_library$sequence == x)
                                             
                                             if (length(idx) > 0) {
                                               indels_library$psi[idx]
                                             } else {
                                               NA
                                             }
                                             
                                           }),
                      library_psi_max = sapply(X = experimental_validations$sequence,
                                               FUN = function(x){
                                                 idx <- which(indels_library$sequence == x)
                                                 
                                                 if (length(idx) > 0) {
                                                   indels_library$psi_max[idx]
                                                 } else {
                                                   NA
                                                 }
                                                 
                                               }),
                      library_psi_min = sapply(X = experimental_validations$sequence,
                                               FUN = function(x){
                                                 idx <- which(indels_library$sequence == x)
                                                 
                                                 if (length(idx) > 0) {
                                                   indels_library$psi_min[idx]
                                                 } else {
                                                   NA
                                                 }
                                                 
                                               }),
                      new_ss = sapply(X =  experimental_validations$SampleID,
                                      FUN = function(x){
                                        x %in% exons_with_new_ss
                                      }),
                      mutation_type = "indel"
                      )



# from 2016 nat comms paper
julien_psi_values <- c("WT" = 49.1,
                       "T21C" = 75.1,
                       "T21G" = 2.7,
                       "C30T" = 23.8,
                       "C36T" = 27.0,
                       "C33T" = 4.6,
                       "C24G" = 5.0,
                       "C24A" = 68.8,
                       "T37G" = 79.9,
                       "T39A" = 53.1,
                       "A50G" = 54.8,
                       "T52C" = 61.8,
                       "C48T" = 21.8,
                       "C49A" = 57.9,
                       "A50C" = 8.2,
                       "C43G" = 70.1,
                       "T59G" = 43.9,
                       "A12C" = 61.1,
                       "G20T" = 19.8,
                       "G20A" = 77.3,
                       "A45G" = 65.3,
                       "A45C" = 7.5,
                       "A44T" = 44.0,
                       "A44G" = 85.1,
                       "C33A" = 5.4,
                       "C33G" = 12.4)


julien_sd_values <- c("WT" = 0.3,
                      "T21C" = 0.5,
                      "T21G" = 0.3,
                      "C30T" = 0.6,
                      "C36T" = 0.3,
                      "C33T" = 0.2,
                      "C24G" = 0.4,
                      "C24A" = 1.5,
                      "T37G" = 1.1,
                      "T39A" = 0.5,
                      "A50G" = 1.5,
                      "T52C" = 1.6,
                      "C48T" = 1.1,
                      "C49A" = 0.7,
                      "A50C" = 0.4,
                      "C43G" = 1.3,
                      "T59G" = 0.8,
                      "A12C" = 1.3,
                      "G20T" = 2,
                      "G20A" = 3.2,
                      "A45G" = 0.3,
                      "A45C" = 1.3,
                      "A44T" = 2,
                      "A44G" = 2.3,
                      "C33A" = 0.1,
                      "C33G" = 1.7)

substitutions_tibble <- tibble(library_id = names(julien_psi_values),
                               experiment_id = names(julien_psi_values),
                               experimental_psi = julien_psi_values,
                               experimental_psi_max = sapply(X = names(julien_psi_values),
                                                             FUN = function(x){
                                                               
                                                               idx <- which(names(julien_psi_values) == x)
                                                               
                                                               psi <- julien_psi_values[idx] + 1.96*(julien_sd_values[idx])
                                                               
                                                               if (psi > 100) {
                                                                 psi <- 100
                                                               }
                                                               
                                                               psi
                                                               
                                                             }),
                               experimental_psi_min = sapply(X = names(julien_psi_values),
                                                             FUN = function(x){
                                                               idx <- which(names(julien_psi_values) == x)
                                                               
                                                               psi <- julien_psi_values[idx] - 1.96*(julien_sd_values[idx])
                                                               
                                                               if (psi < 0) {
                                                                 psi <- 0
                                                               }
                                                               
                                                               psi
                                                               
                                                             }),
                               
                               library_psi = sapply(X = names(julien_psi_values),
                                                    FUN = function(x){
                                                      idx <- which(indels_library$id == x)
                                                      
                                                      psi <- indels_library$psi[idx]
                                                      
                                                      psi
                                                      
                                                    }),
                               library_psi_max = sapply(X = names(julien_psi_values),
                                                        FUN = function(x){
                                                          idx <- which(indels_library$id == x)
                                                          
                                                          psi <- indels_library$psi_max[idx]
                                                          
                                                          psi
                                                          
                                                        }),
                               
                               library_psi_min = sapply(X = names(julien_psi_values),
                                                        FUN = function(x){
                                                          idx <- which(indels_library$id == x)
                                                          
                                                          psi <- indels_library$psi_min[idx]
                                                          
                                                          psi
                                                          
                                                        }),
                               
                               new_ss = NA,
                               
                               mutation_type = "substitution"
                               
                               )





indels_tibble$experimental_psi <- sapply(X = indels_tibble$experimental_psi,
                                         FUN = function(x){
                                           
                                           starting_psi <- (75.7 + 68.4)/2
                                           delta_psi <- x - starting_psi
                                           
                                           A <- (starting_psi^2 - 100*starting_psi + delta_psi*starting_psi - 100*delta_psi) / (starting_psi*(delta_psi + starting_psi - 100))
                                           
                                           new_starting_psi <- 49.1
                                           
                                           new_delta_psi <- (  (100 * A * new_starting_psi) / (100 - new_starting_psi + A*new_starting_psi)  ) - new_starting_psi
                                           
                                           new_psi <- new_delta_psi + new_starting_psi
                                           
                                           new_psi
                                           
                                         }) 


indels_tibble$experimental_psi_max <- sapply(X = indels_tibble$experimental_psi_max,
                                             FUN = function(x){
                                               
                                               starting_psi <- (75.7 + 68.4)/2
                                               delta_psi <- x - starting_psi
                                               
                                               A <- (starting_psi^2 - 100*starting_psi + delta_psi*starting_psi - 100*delta_psi) / (starting_psi*(delta_psi + starting_psi - 100))
                                               
                                               new_starting_psi <- 49.1
                                               
                                               new_delta_psi <- (  (100 * A * new_starting_psi) / (100 - new_starting_psi + A*new_starting_psi)  ) - new_starting_psi
                                               
                                               new_psi <- new_delta_psi + new_starting_psi
                                               
                                               new_psi
                                               
                                             }) 


indels_tibble$experimental_psi_min <- sapply(X = indels_tibble$experimental_psi_min,
                                             FUN = function(x){
                                               
                                               starting_psi <- (75.7 + 68.4)/2
                                               delta_psi <- x - starting_psi
                                               
                                               A <- (starting_psi^2 - 100*starting_psi + delta_psi*starting_psi - 100*delta_psi) / (starting_psi*(delta_psi + starting_psi - 100))
                                               
                                               new_starting_psi <- 49.1
                                               
                                               new_delta_psi <- (  (100 * A * new_starting_psi) / (100 - new_starting_psi + A*new_starting_psi)  ) - new_starting_psi
                                               
                                               new_psi <- new_delta_psi + new_starting_psi
                                               
                                               new_psi
                                               
                                             }) 




plot_tibble <- rbind(indels_tibble,
                     substitutions_tibble)


ggplot(data = plot_tibble %>%
         filter(experiment_id != "WT"),
       mapping = aes(y = experimental_psi,
                     x = library_psi)) +
  geom_smooth(method = "lm",
              fill = "#F9C45F",
              colour = "#F9C45F",
              alpha = 0.2) +
  geom_errorbar(mapping = aes(ymin = experimental_psi_min,
                              ymax = experimental_psi_max),
                width = 0,
                position=position_dodge(.9),
                colour = "gray30") +
  geom_errorbar(mapping = aes(xmin = library_psi_min,
                              xmax = library_psi_max),
                width = 0,
                position=position_dodge(.9),
                colour = "gray30") +
  geom_point(size = 3,
             colour = "gray30",
             mapping = aes(shape = mutation_type)) +
  coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  ylab("PSI (experimental)") +
  xlab("PSI (RNAseq)")  +
  scale_shape_manual(name = "mutation type", values = c(16,15)) +
  ggtitle("rho = 0.91")
  

ggsave(filename = "129_experimental_validation.pdf",
       width = 7,
       height = 7,
       useDingbats = F)


cor.test(x = (plot_tibble %>%
                filter(experiment_id != "WT"))$library_psi,
         y = (plot_tibble %>%
                filter(experiment_id != "WT"))$experimental_psi,
         method = "spearman")


cor.test(x = (plot_tibble %>%
                filter(experiment_id != "WT" & mutation_type == "indel"))$library_psi,
         y = (plot_tibble %>%
                filter(experiment_id != "WT" & mutation_type == "indel"))$experimental_psi,
         method = "spearman")


cor.test(x = (plot_tibble %>%
                filter(experiment_id != "WT" & mutation_type != "indel"))$library_psi,
         y = (plot_tibble %>%
                filter(experiment_id != "WT" & mutation_type != "indel"))$experimental_psi,
         method = "spearman")
