library(tidyverse)
library(scales)


k <- 3

possible_kmers <- sort(apply(X = expand.grid(rep(list(c("A","C","G","T")), k)),
                             MARGIN = 1,
                             FUN = function(x){
                               paste(x, sep = "", collapse = "")
                             }))









# reverse complement function
rev.comp<-function(x,rev=TRUE)
{
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"    
    if(xx[bbb]=="C") y[bbb]<-"G"    
    if(xx[bbb]=="G") y[bbb]<-"C"    
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)  
}


# load table
bed_table <- read.table(file = "../Revisions_eLife/drive-download-20200926T084945Z-001/238_GTEx_Exon_Regions_With_Sequences.bed",
                        header = F,
                        sep = "\t",
                        col.names = c("chromosome",
                                      "start",
                                      "end",
                                      "id",
                                      "strand",
                                      "gene",
                                      "length",
                                      "sequence"),
                        stringsAsFactors = F)


# reverse complement any sequence in the reverse strand
bed_table$sequence <- apply(X = bed_table[,c("sequence", "strand")],
                            MARGIN = 1,
                            FUN = function(x){
                              mySequence <- x[1]
                              myStrand <- x[2]
                              
                              if (myStrand == "-") {
                                mySequence <- rev.comp(mySequence)
                              }
                              
                              mySequence
                            })





kmer_counts_matrix <- t(sapply(X = as.character(bed_table$sequence),
                               FUN = function(x){
                                 sapply(X = possible_kmers,
                                        FUN = function(y){
                                          string <- x
                                          pattern <- y
                                          nstring <- nchar(string)
                                          npattern <- nchar(pattern)
                                          stringlist <- substring(string, 1:(nstring-npattern+1), npattern:nstring)
                                          sum(stringlist==pattern)
                                        })
                               }))

rownames(kmer_counts_matrix) <- as.character(bed_table$id)


















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


All.Mean.PSIs <- c()
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
  
  All.Mean.PSIs <- c(All.Mean.PSIs, unname(Mean.PSI.Values))
  All.Exon.IDs <- c(All.Exon.IDs, names(Mean.PSI.Values))
  All.Tissues <- c(All.Tissues, rep(tissue.no.spaces, length(Mean.PSI.Values)))
}


psi_tibble <- tibble(psi = All.Mean.PSIs,
                     id = All.Exon.IDs,
                     tissue = All.Tissues)

psi_tibble_wide <- psi_tibble %>%
  spread(key = tissue, value = psi)

psi_matrix <- as.matrix(x = psi_tibble_wide %>%
                          select(AdiposeTissue:Vagina))
rownames(psi_matrix) <- as.character(psi_tibble_wide$id)

psi_matrix <- psi_matrix[rownames(kmer_counts_matrix),]













median_psi_per_kmer <- sapply(X = possible_kmers,
                              FUN = function(x){
                                idx <- which(kmer_counts_matrix[,x]>9)
                                median(psi_matrix[idx,"AdiposeTissue"], na.rm = T)
                              })




load("109_median_scores.RData")





gtex_order <- sapply(X = possible_kmers,
                     FUN = function(x){
                       gtex <- possible_kmers[order(median_psi_per_kmer, decreasing = T)]
                       which(gtex == x)
                     })

library_order <- sapply(X = possible_kmers,
                        FUN = function(x){
                          library <- names(median_scores_library)
                          which(library == x)
                        })



gtex_order_3mer <- gtex_order/length(gtex_order)
library_order_3mer <- library_order/length(library_order)
possible_3mers <- possible_kmers

correlation_3mers <- cor.test(x = library_order, y = gtex_order, method = "spearman")




k <- 2

possible_kmers <- sort(apply(X = expand.grid(rep(list(c("A","C","G","T")), k)),
                             MARGIN = 1,
                             FUN = function(x){
                               paste(x, sep = "", collapse = "")
                             }))




kmer_counts_matrix <- t(sapply(X = as.character(bed_table$sequence),
                               FUN = function(x){
                                 sapply(X = possible_kmers,
                                        FUN = function(y){
                                          string <- x
                                          pattern <- y
                                          nstring <- nchar(string)
                                          npattern <- nchar(pattern)
                                          stringlist <- substring(string, 1:(nstring-npattern+1), npattern:nstring)
                                          sum(stringlist==pattern)
                                        })
                               }))

rownames(kmer_counts_matrix) <- as.character(bed_table$id)





median_psi_per_kmer <- sapply(X = possible_kmers,
                              FUN = function(x){
                                idx <- which(kmer_counts_matrix[,x]>19)
                                median(psi_matrix[idx,"AdiposeTissue"], na.rm = T)
                              })

load("110_median_scores.RData")





gtex_order <- sapply(X = possible_kmers,
                     FUN = function(x){
                       gtex <- possible_kmers[order(median_psi_per_kmer, decreasing = T)]
                       which(gtex == x)
                     })

library_order <- sapply(X = possible_kmers,
                        FUN = function(x){
                          library <- names(median_scores_library)
                          which(library == x)
                        })



gtex_order_2mer <- gtex_order/length(gtex_order)
library_order_2mer <- library_order/length(library_order)
possible_2mers <- possible_kmers

correlation_2mers <- cor.test(x = library_order, y = gtex_order, method = "spearman")









####

scatterplot_tibble <- tibble(gtex = c(gtex_order_3mer, gtex_order_2mer),
                             library = c(library_order_3mer, library_order_2mer),
                             kmer = c(possible_3mers, possible_2mers),
                             kmer_size = c(rep(3,length(gtex_order_3mer)),
                                           rep(2, length(gtex_order_2mer))))
scatterplot_tibble$kmer <- gsub(pattern = "T",
                                replacement = "U",
                                x = scatterplot_tibble$kmer)



ggplot(data = scatterplot_tibble,
       mapping = aes(x = 1-gtex,
                     y = 1-library,
                     label = kmer,
                     colour = as.factor(kmer_size))) +
  geom_abline(intercept = 0, slope = 1, colour = "gray60", lty = 2, size = 1) +
  geom_text() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  scale_colour_manual(values = c("#F9C45F", "black")) +
  xlab("normalised rank (GTEx adipose tissue)") +
  ylab("normalised rank (mutation library)")



ggsave(filename = "123_kmer_rank_correlations_adipose_tissue.pdf",
       height = 6,
       width = 6,
       useDingbats = FALSE)

###








