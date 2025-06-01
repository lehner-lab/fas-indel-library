
library(ggplot2)
source("../Revisions_eLife/gtex/185_Functions.R")

Tissues <- c("Adipose Tissue",
             "Adrenal Gland",
             #"Bladder",
             "Blood",
             "Blood Vessel",
             #"Bone Marrow",
             "Brain",
             "Breast",
             #"Cervix Uteri",
             "Colon",
             "Esophagus",
             #"Fallopian Tube",
             "Heart",
             "Kidney",
             "Liver",
             "Lung",
             "Muscle",
             "Nerve",
             "Ovary",
             "Pancreas",
             "Pituitary",
             "Prostate",
             "Salivary Gland",
             "Skin",
             "Small Intestine",
             "Spleen",
             "Stomach",
             "Testis",
             "Thyroid",
             "Uterus",
             "Vagina")


All.Starting.PSIs <- c()
All.Exon.IDs <- c()
All.Tissues <- c()

for (each.tissue in Tissues) {
  print(each.tissue)
  
  tissue.no.spaces <- gsub(pattern = " ",
                           replacement = "",
                           x = each.tissue,
                           fixed = TRUE)
  
  file.location <- paste("../Revisions_eLife/gtex/182_GTEx_PSI_Distributions_Plus_IDs/ALL_PSIandID_",
                         tissue.no.spaces,
                         ".RData",
                         sep = "")
  
  load(file.location)
  
  
  
  if (length(Mean.PSI.Values)==0){
    next
  }
  
  All.Starting.PSIs <- c(All.Starting.PSIs, unname(Mean.PSI.Values))
  All.Exon.IDs <- c(All.Exon.IDs, names(Mean.PSI.Values))
  All.Tissues <- c(All.Tissues, rep(tissue.no.spaces, length(Mean.PSI.Values)))
}



Plot.DF <- data.frame(Starting.PSI = All.Starting.PSIs,
                      Tissue = All.Tissues,
                      Exon.ID = All.Exon.IDs)


Plot.DF$length <- sapply(X = Plot.DF$Exon.ID,
                         FUN = function(x){
                           abs(as.numeric(strsplit(x, "_")[[1]][5]) - as.numeric(strsplit(x, "_")[[1]][6])) + 1
                         })


Plot.DF$length_group <- findInterval(x = Plot.DF$length, vec = seq(0,100,10), all.inside = T)

ggplot(data = Plot.DF[which(Plot.DF$Tissue == "Heart"),],
       mapping = aes(x = length_group,
                     y = Starting.PSI,
                     group = length_group)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid = element_blank())


unique(Plot.DF$length_group)
library(tidyverse)

table(Plot.DF$Tissue)

proportions_matrix <- sapply(X = 1:10,
                             FUN = function(x){
                               
                               sapply(X = seq(0,0.8,0.2),
                                      FUN = function(y){
                                        if (y == 0.8) {
                                          mean((Plot.DF$Starting.PSI[ Plot.DF$length_group == x & Plot.DF$Tissue == "AdiposeTissue"] <= y+0.2) & (Plot.DF$Starting.PSI[ Plot.DF$length_group == x  & Plot.DF$Tissue == "AdiposeTissue"] >= y) , na.rm = T)
                                        } else {
                                          mean((Plot.DF$Starting.PSI[ Plot.DF$length_group == x & Plot.DF$Tissue == "AdiposeTissue"] < y+0.2) & (Plot.DF$Starting.PSI[ Plot.DF$length_group == x  & Plot.DF$Tissue == "AdiposeTissue"] >= y) , na.rm = T)
                                        }
                                      })
                             })

proportions_vector <- c()

for (i in 1:ncol(proportions_matrix)) {
  proportions_vector <- c(proportions_vector,
                          proportions_matrix[,i])
}

plot_tibble <- tibble(length_group = rep(1:10, each = 5),
                      inclusion_group = rep(1:5, 10),
                      proportion = proportions_vector)




colour_function <- colorRampPalette(c("white","#F9C45F"))
my_colours <- colour_function(5)


ggplot(data = plot_tibble,
       mapping = aes(x = factor(length_group),
                     y = proportion,
                     fill = factor(inclusion_group, levels = 10:1))) +
  geom_bar(stat = "identity", position = "stack", colour = "black") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1.5) +
  scale_fill_manual(values = my_colours)


ggsave(filename = "117_gtex_psi_vs_exon_length_adipose_tissue.pdf",
       height = 10,
       width = 10)



write_delim(x = plot_tibble,
            delim = "\t",
            file = "117_figure_03b.tsv")
