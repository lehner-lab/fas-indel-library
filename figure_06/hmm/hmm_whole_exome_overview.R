
setwd("/omics/groups/OE0433/internal/pablo/projects/sandbox/illumina/")

library(tidyverse)




load("hmm_hidden_states_split_by_exon.RData")
load("hmm_exon_lengths.RData")






load("001_master_bed_table_exons_in_80pcnt_samples_with_refseqs_and_psiVals.RData")




# only look at exons of length 50-150
# which(grepl(pattern = "FAS", x = bed_file_exons_in_80_pcnt_of_samples$exon_id) & bed_file_exons_in_80_pcnt_of_samples$length == 63)
# 
# idx <- c(11745, which(bed_file_exons_in_80_pcnt_of_samples$length == 100))
idx <- which(bed_file_exons_in_80_pcnt_of_samples$length >= 50 & bed_file_exons_in_80_pcnt_of_samples$length <= 200)
bed_file_exons_in_80_pcnt_of_samples <- bed_file_exons_in_80_pcnt_of_samples[idx,]

my_exons <- bed_file_exons_in_80_pcnt_of_samples$exon_id






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

    next
  }
  
  # will now get the mean psi for this exon
  bed_idx <- which(bed_file_exons_in_80_pcnt_of_samples$exon_id == this_exon)
  this_exon_psi <- 100 * mean(x = suppressWarnings(as.numeric(bed_file_exons_in_80_pcnt_of_samples[bed_idx,9:ncol(bed_file_exons_in_80_pcnt_of_samples)])), na.rm = TRUE)
  wt_psi_vector <- c(wt_psi_vector,
                     this_exon_psi)
  
  
}



names(wt_psi_vector) <- my_exons
wt_psi_vector <- wt_psi_vector[-which(is.na(wt_psi_vector))]











tibble_exome_overview <- tibble(psi = wt_psi_vector,
                                  percent_E = sapply(X = hidden_states_split_by_exon,
                                                     FUN = function(x){
                                                       100 * mean(x == 1)
                                                     }),
                                  sum_E = sapply(X = hidden_states_split_by_exon,
                                                 FUN = function(x){
                                                   sum(x == 1)
                                                 }),
                                  
                                  percent_N = sapply(X = hidden_states_split_by_exon,
                                                     FUN = function(x){
                                                       100 * mean(x == 2)
                                                     }),
                                  sum_N = sapply(X = hidden_states_split_by_exon,
                                                 FUN = function(x){
                                                   sum(x == 2)
                                                 }),
                                  
                                  percent_S = sapply(X = hidden_states_split_by_exon,
                                                     FUN = function(x){
                                                       100 * mean(x == 3)
                                                     }),
                                  sum_S = sapply(X = hidden_states_split_by_exon,
                                                 FUN = function(x){
                                                   sum(x == 3)
                                                 }),
                                  length = exon_lengths)

tibble_exome_overview$exon = names(hidden_states_split_by_exon)


colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 7, palette = "Spectral")[2:5]))
my_colours <- colour_function(11)


# 
# ggplot(data = tibble_exome_overview,
#        mapping = aes(x = percent_S,
#                      y = sum_S,
#                      colour = psi)) +
#   geom_point(size = 1) +
#   scale_colour_gradientn(colours = my_colours, 
#                          values = scales::rescale(c(0,
#                                                     10,
#                                                     20,
#                                                     30,
#                                                     40,
#                                                     50,
#                                                     60,
#                                                     70,
#                                                     80,
#                                                     90,
#                                                     100)),
#                          limits = c(0,100),
#                          guide = "colourbar",
#                          na.value = "black",
#                          oob = scales::squish)
# theme_bw() +
#   theme(panel.grid = element_blank())



ggplot(data = tibble_exome_overview,
       mapping = aes(x = percent_S,
                     y = percent_E,
                     colour = psi)) +
  geom_point(size = 2) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(c(0,100)) +
  ylim(c(0,100))


ggtern::ggtern(data = tibble_exome_overview,
               mapping = aes(x = percent_S,
                             y = percent_E,
                             z = percent_N,
                             colour = psi)) +
  geom_point(size = 2) +
  ggtern::theme_rgbw() +
  theme_bw() +
  ggtern::theme_gridsontop() +
  theme(panel.grid.major = element_line(color = "black"),
        panel.grid.minor = element_line(color = "black")) +
  scale_color_gradient2(midpoint = 50, low = "#42B9C4", high = "#F9C45F") +
  # theme_bw() +
  # theme(panel.grid = element_blank()) +
  xlim(c(0,100)) +
  ylim(c(0,100)) +
  xlab("% S state") +
  ylab("% E state") +
  ggtern::zlab("% N state")

# exported manually to pdf 7x7


ggplot(data = tibble_exome_overview,
       mapping = aes(x = percent_S,
                     y = percent_E,
                     colour = psi)) +
  geom_point(size = 2) +
  scale_color_gradient2(midpoint = 50, low = "#42B9C4", high = "#F9C45F") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(c(0,100)) +
  ylim(c(0,100))


# 
# 
# vector_percent_s <- c()
# vector_percent_e <- c()
# 
# 
# for (i in seq(0,100, 0.5)) {
#   
#   for (j in seq(0,100,0.5)) {
#     
#     if (i+j > 100) {
#       next
#     }
#     
#     vector_percent_s <- c(vector_percent_s,
#                           i)
#     vector_percent_e <- c(vector_percent_e,
#                           j)
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# 
# loess_model <- loess(psi ~ percent_S + percent_E,
#                      data = tibble_exome_overview,
#                      span = 0.75, family = "symmetric")
# 
# 
# model_predictions <- predict(object = loess_model,
#                              newdata = tibble(percent_S = vector_percent_s,
#                                               percent_E = vector_percent_e))
# 
# 
# 
# plot_tibble <- tibble(percent_s = vector_percent_s,
#                       percent_e = vector_percent_e,
#                       psi = model_predictions)
# 
# library(metR)
# 
# ggplot(data = plot_tibble,
#        mapping = aes(x = percent_s,
#                      y = percent_e,
#                      z = psi)) +
#   geom_contour_fill(na.fill = T)+
#   geom_contour(color = "black") +
#   scale_fill_gradient2(midpoint = 50, low = "#42B9C4", high = "#F9C45F", name = "PSI") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         aspect.ratio = 1) +
#   geom_abline(slope = -1,
#               intercept = 100) +
#   xlab("% S state") +
#   ylab("% E state") +
#   scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,100))
# 
# 
# 
# 
# 
# 
# 
# vector_percent_s <- c()
# vector_percent_e <- c()
# vector_percent_n <- c()
# 
# 
# for (i in seq(0,100, 1)) {
#   
#   for (j in seq(0,100,1)) {
#     
#     for (k in seq(0,100,1)) {
#       
#       if (i+j+k > 100) {
#         next
#       }
#       vector_percent_s <- c(vector_percent_s,
#                             i)
#       vector_percent_e <- c(vector_percent_e,
#                             j)
#       
#       vector_percent_n <- c(vector_percent_n,
#                             k)
#     }
#    
#     
#    
#     
#   }
#   
# }
# 
# 
# loess_model <- loess(psi ~ percent_S + percent_E + percent_N,
#                      data = tibble_exome_overview,
#                      span = 0.75, family = "symmetric")
# 
# 
# model_predictions <- predict(object = loess_model,
#                              newdata = tibble(percent_S = vector_percent_s,
#                                               percent_E = vector_percent_e,
#                                               percent_N = vector_percent_n))
# 
# 
# 
# plot_tibble <- tibble(percent_s = vector_percent_s,
#                       percent_e = vector_percent_e,
#                       percent_n = vector_percent_n,
#                       psi = model_predictions)
# 
# 
# 
# 
# ggtern::ggtern(data = plot_tibble,
#                mapping = aes(x = percent_s,
#                              y = percent_e,
#                              z = percent_n,
#                              colour = psi)) +
#   geom_point(size = 2) +
#   ggtern::theme_rgbw() +
#   theme_bw() +
#   ggtern::theme_gridsontop() +
#   theme(panel.grid.major = element_line(color = "black"),
#         panel.grid.minor = element_line(color = "black")) +
#   scale_color_gradient2(midpoint = 50, low = "#42B9C4", high = "#F9C45F") +
#   # theme_bw() +
#   # theme(panel.grid = element_blank()) +
#   xlim(c(0,100)) +
#   ylim(c(0,100)) +
#   xlab("% S state") +
#   ylab("% E state") +
#   ggtern::zlab("% N state")











tibble_exome_overview$length2 <- sapply(X = tibble_exome_overview$length,
                                        FUN = function(x){
                                          if (x > 200) {
                                            x <- 200
                                          }
                                          x
                                        })
ggplot(data = tibble_exome_overview,
       mapping = aes(x = percent_E,
                     y = percent_S,
                     colour = length2)) +
  geom_point(size = 2) +
  scale_color_gradient2(midpoint = 100, low = "#42B9C4", high = "#F9C45F") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlim(c(0,100)) +
  ylim(c(0,100))




ggtern::ggtern(data = tibble_exome_overview,
               mapping = aes(x = percent_S,
                             y = percent_E,
                             z = percent_N,
                             colour = length2)) +
  geom_point(size = 2) +
  ggtern::theme_rgbw() +
  theme_bw() +
  ggtern::theme_gridsontop() +
  theme(panel.grid.major = element_line(color = "black"),
        panel.grid.minor = element_line(color = "black")) +
  scale_color_gradient2(midpoint = 110, low = "#42B9C4", high = "#F9C45F") +
  # theme_bw() +
  # theme(panel.grid = element_blank()) +
  xlim(c(0,100)) +
  ylim(c(0,100)) +
  xlab("% S state") +
  ylab("% E state") +
  ggtern::zlab("% N state")

# gives an error but I can manually export to pdf
ggsave(filename = "hmm_exome_length_tern.pdf",
       width = 7,
       height = 7,
       useDingbats = F)




# loess_model <- loess(length2 ~ percent_S + percent_E,
#                      data = tibble_exome_overview,
#                      span = 0.75, family = "symmetric")
# 
# 
# model_predictions <- predict(object = loess_model,
#                              newdata = tibble(percent_S = vector_percent_s,
#                                               percent_E = vector_percent_e))
# 
# 
# 
# plot_tibble <- tibble(percent_s = vector_percent_s,
#                       percent_e = vector_percent_e,
#                       length = model_predictions)
# 
# library(metR)
# 
# 
# ggplot(data = plot_tibble,
#        mapping = aes(x = percent_s,
#                      y = percent_e,
#                      z = length)) +
#   geom_contour_fill(na.fill = T)+
#   geom_contour(color = "black") +
#   scale_fill_gradient2(midpoint = 100, low = "#42B9C4", high = "#F9C45F", name = "length") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         aspect.ratio = 1) +
#   geom_abline(slope = -1,
#               intercept = 100) +
#   xlab("% S state") +
#   ylab("% E state") +
#   scale_x_continuous(expand = c(0,0), limits = c(0,100)) +
#   scale_y_continuous(expand = c(0,0), limits = c(0,100))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ggplot(data = tibble_exome_overview,
#        mapping = aes(x = percent_S,
#                      y = percent_E,
#                      colour = sum_S + sum_E)) +
#   geom_point(size = 2) +
#   scale_color_gradient2(midpoint = 100, low = "#42B9C4", high = "#F9C45F") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlim(c(0,100)) +
#   ylim(c(0,100))
# 
# 
# ggplot(data = tibble_exome_overview,
#        mapping = aes(x = percent_E,
#                      y = percent_S,
#                      colour = findInterval(x = length2, vec = c(0,100,150,200), all.inside = T))) +
#   geom_density_2d(h = c(50,50), inherit.aes = T) +
#   # geom_point(size = 2) +
#   # scale_color_gradient2(midpoint = 100, low = "#42B9C4", high = "#F9C45F") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlim(c(0,100)) +
#   ylim(c(0,100))
