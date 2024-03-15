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






####



scatterplot_tibble <- tibble(gtex = gtex_order_3mer,
                             library = library_order_3mer,
                             kmer = possible_3mers,
                             kmer_size = rep(3,length(gtex_order_3mer)))
scatterplot_tibble$kmer <- gsub(pattern = "T",
                                replacement = "U",
                                x = scatterplot_tibble$kmer)


proline <- c("CCU", "CCC", "CCA", "CCG") # Pro
# glycine <- c("GGU", "GGC", "GGA", "GGG") # Gly

nonpolar_short <- c("GGU", "GGC", "GGA", "GGG", # Gly
                    "GCU", "GCC", "GCA", "GCG", # Ala
                    "GUU", "GUC", "GUA", "GUG" # Val
                    )
nonpolar <- c("UUA", "UUG", "CUU", "CUC", "CUA", "CUG", # Leu
              "AUU", "AUC", "AUA", #Ile
              "AUG" #Met
              )

aromatic <- c("UUU", "UUC", # Phe
              "UGG", # Trp
              "UAU", "UAC" # Tyr
              
)

polar <- c("UCU", "UCC", "UCA", "UCG", # Ser
           "ACU", "ACC", "ACA", "ACG", # Thr
           
           "CAA", "CAG", # Gln
           "AAU", "AAC", # Asn
           "UGU", "UGC", # Cys
           "AGU", "AGC") # Ser

basic <- c("CAU", "CAC", # His
           "AAA", "AAG", # Lys
           "CGU", "CGC", "CGA", "CGG", # Arg
           "AGA", "AGG") # Arg

acidic <- c("GAU", "GAC", # Asp
            "GAA", "GAG") # Glu

stops <- c("UAA", "UAG", "UGA")



charged <- unique(c(basic, acidic))
nonpolar_all <- unique(c(nonpolar_short, nonpolar, proline))





# better combine all nonpolar together and all charged together
prop.table(x = table(scatterplot_tibble$codon_class, scatterplot_tibble$uc_counts),
           margin = 2)[c("nonpolar", "aromatic", "polar", "charged", "stop"),]

pdf(file = "123b_barplot_aminoacid_vs_uc_content.pdf",
    width = 8,
    height = 8,
    useDingbats = F)
barplot(prop.table(x = table(scatterplot_tibble$codon_class, scatterplot_tibble$uc_counts),
                   margin = 2)[c("nonpolar", "aromatic", "polar", "charged", "stop"),]*100,
        las = 1,
        border = NA,
        main = "amino acid class vs UC content")
dev.off()




genetic_code <- tibble(codon = c("UUU", "UUC", "UUA", "UUG",
                                 "CUU", "CUC", "CUA", "CUG",
                                 "AUU", "AUC", "AUA", "AUG",
                                 "GUU", "GUC", "GUA", "GUG",
                                 "UCU", "UCC", "UCA", "UCG",
                                 "CCU", "CCC", "CCA", "CCG",
                                 "ACU", "ACC", "ACA", "ACG",
                                 "GCU", "GCC", "GCA", "GCG",
                                 "UAU", "UAC", "UAA", "UAG",
                                 "CAU", "CAC", "CAA", "CAG",
                                 "AAU", "AAC", "AAA", "AAG",
                                 "GAU", "GAC", "GAA", "GAG",
                                 "UGU", "UGC", "UGA", "UGG",
                                 "CGU", "CGC", "CGA", "CGG",
                                 "AGU", "AGC", "AGA", "AGG",
                                 "GGU", "GGC", "GGA", "GGG"),
                       
                       aa3 = c("Phe", "Phe", "Leu", "Leu",
                               "Leu", "Leu", "Leu", "Leu",
                               "Ile", "Ile", "Ile", "Met",
                               "Val", "Val", "Val", "Val",
                               "Ser", "Ser", "Ser", "Ser",
                               "Pro", "Pro", "Pro", "Pro",
                               "Thr", "Thr", "Thr", "Thr",
                               "Ala", "Ala", "Ala", "Ala",
                               "Tyr", "Tyr", "Stop","Stop",
                               "His", "His", "Gln", "Gln",
                               "Asn", "Asn", "Lys", "Lys",
                               "Asp", "Asp", "Glu", "Glu",
                               "Cys", "Cys", "Stop","Trp",
                               "Arg", "Arg", "Arg", "Arg",
                               "Ser", "Ser", "Arg", "Arg",
                               "Gly", "Gly", "Gly", "Gly"),
                       
                       aa1 = c("F", "F", "L", "L",
                               "L", "L", "L", "L",
                               "I", "I", "I", "M",
                               "V", "V", "V", "V",
                               "S", "S", "S", "S",
                               "P", "P", "P", "P",
                               "T", "T", "T", "T",
                               "A", "A", "A", "A",
                               "Y", "Y", "*", "*",
                               "H", "H", "Q", "Q",
                               "N", "N", "K", "K",
                               "D", "D", "E", "E",
                               "C", "C", "*", "W",
                               "R", "R", "R", "R",
                               "S", "S", "R", "R",
                               "G", "G", "G", "G"))





?Peptides::hydrophobicity
Peptides::hydrophobicity("H", scale = "KyteDoolittle")


scatterplot_tibble$codon_class <- sapply(X = scatterplot_tibble$kmer,
                                         FUN = function(x){
                                           result <- NA
                                           if (x %in% nonpolar_all) {
                                             result <- "nonpolar"
                                           }
                                           if (x %in% polar) {
                                             result <- "polar"
                                           }
                                           if (x %in% charged) {
                                             result <- "charged"
                                           }
                                           if (x %in% stops) {
                                             result <- "stop"
                                           }
                                           if (x %in% aromatic) {
                                             result <- "aromatic"
                                           }
                                           
                                           result
                                         })

scatterplot_tibble$codon_class_2 <- sapply(X = scatterplot_tibble$kmer,
                                         FUN = function(x){
                                           result <- NA
                                           if (x %in% nonpolar) {
                                             result <- "nonpolar long"
                                           }
                                           if (x %in% nonpolar_short) {
                                             result <- "nonpolar short"
                                           }
                                           if (x %in% proline) {
                                             result <- "proline"
                                           }
                                           if (x %in% polar) {
                                             result <- "polar"
                                           }
                                           if (x %in% charged) {
                                             result <- "charged"
                                           }
                                           if (x %in% stops) {
                                             result <- "stop"
                                           }
                                           if (x %in% aromatic) {
                                             result <- "aromatic"
                                           }
                                           
                                           result
                                         })


scatterplot_tibble$hydropathy <- sapply(X = scatterplot_tibble$kmer,
                                         FUN = function(x){
                                           this_idx <- which(genetic_code$codon == x)
                                           this_aa1 <- genetic_code$aa1[this_idx]
                                           if (this_aa1 == "*") {
                                             NA
                                           } else {
                                             Peptides::hydrophobicity(this_aa1, scale = "KyteDoolittle")
                                           }
                                         })



ggplot(data = scatterplot_tibble,
       mapping = aes(x = 1-gtex,
                     y = 1-library,
                     label = kmer,
                     colour = hydropathy)) +
  geom_abline(intercept = 0, slope = 1, colour = "gray60", lty = 2, size = 1) +
  geom_text() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  # scale_colour_manual(values = c("#F9C45F", "black")) +
  xlab("normalised rank (GTEx adipose tissue)") +
  ylab("normalised rank (mutation library)")




ggplot(data = scatterplot_tibble,
       mapping = aes(x = 1-gtex,
                     y = 1-library,
                     label = kmer,
                     colour = as.factor(codon_class_2))) +
  geom_abline(intercept = 0, slope = 1, colour = "gray60", lty = 2, size = 1) +
  geom_text() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  # scale_colour_manual(values = c("#F9C45F", "black")) +
  xlab("normalised rank (GTEx adipose tissue)") +
  ylab("normalised rank (mutation library)")


ggplot(data = scatterplot_tibble,
       mapping = aes(x = 1-library,
                     label = kmer,
                     colour = as.factor(codon_class))) +
  geom_density() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  # scale_colour_manual(values = c("#F9C45F", "black")) +
  xlab("normalised rank (GTEx adipose tissue)") +
  ylab("normalised rank (mutation library)")




ggplot(data = scatterplot_tibble,
       mapping = aes(x = 1-gtex,
                     label = kmer,
                     colour = as.factor(codon_class))) +
  geom_density() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  # scale_colour_manual(values = c("#F9C45F", "black")) +
  xlab("normalised rank (GTEx adipose tissue)") +
  ylab("normalised rank (mutation library)")





ggsave(filename = "123b_kmer_rank_correlations_adipose_tissue.pdf",
       height = 6,
       width = 6,
       useDingbats = FALSE)

###

str_count()

str_count(q.data$string, "a")

scatterplot_tibble$uc_counts <- sapply(X = scatterplot_tibble$kmer,
                                       FUN = function(x){
                                         str_count(string = x, pattern = "U") + str_count(string = x, pattern = "C")
                                       })


for (i in 0:3) {
  
}


