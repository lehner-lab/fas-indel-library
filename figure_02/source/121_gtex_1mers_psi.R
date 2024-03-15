library(tidyverse)
library(scales)


k <- 1

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
                                          sum(stringlist==pattern) / nchar(string)
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




for (this_nucleotide in possible_kmers) {
  proportions_matrix <- sapply(X = seq(0,0.5,0.1), # % of kmer
                               FUN = function(x){
                                 
                                 if (x == 0.5) {
                                   idx <- which(kmer_counts_matrix[,this_nucleotide]<= 1 & kmer_counts_matrix[,this_nucleotide]>= x)
                                 } else {
                                   idx <- which(kmer_counts_matrix[,this_nucleotide]< x+0.1 & kmer_counts_matrix[,this_nucleotide]>= x)
                                 }
                                 
                                 sapply(X = seq(0,0.8,0.2), # PSI
                                        FUN = function(y){
                                          if (y == 0.8) {
                                            mean((psi_matrix[idx,"AdiposeTissue"] <= y+0.2) & (psi_matrix[idx,"AdiposeTissue"] >= y) , na.rm = T)
                                          } else {
                                            mean((psi_matrix[idx,"AdiposeTissue"] < y+0.2) & (psi_matrix[idx,"AdiposeTissue"] >= y) , na.rm = T)
                                          }
                                        })
                               })
  
  proportions_vector <- c()
  
  for (i in 1:ncol(proportions_matrix)) {
    proportions_vector <- c(proportions_vector,
                            proportions_matrix[,i])
  }
  
  plot_tibble <- tibble(percentage_group = rep(1:ncol(proportions_matrix), each = 5),
                        inclusion_group = rep(1:5, ncol(proportions_matrix)),
                        proportion = proportions_vector)
  
  
  colour_function <- colorRampPalette(c("white","#F9C45F"))
  my_colours <- colour_function(5)
  
  
  ggplot(data = plot_tibble,
         mapping = aes(x = percentage_group,
                       y = proportion,
                       fill = factor(inclusion_group, levels = 5:1))) +
    geom_bar(stat = "identity", position = "stack", colour = "black") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.ticks.x=element_blank()
    ) +
    scale_fill_manual(values = my_colours,
                      name = "inclusion") +
    xlab("") +
    scale_x_continuous(breaks = 1:6,
                       labels=c("0–9",
                                "10–19",
                                "20–29",
                                "30–39",
                                "40–49",
                                "50+")) +
    ggtitle(this_nucleotide)
  
  file_name <- paste("121_gtex_adipose_tissue_",
                     this_nucleotide,
                     ".pdf",
                     sep = "",
                     collapse = "")
  ggsave(filename = file_name,
         height = 10,
         width = 10,
         useDingbats = FALSE)
  
}




