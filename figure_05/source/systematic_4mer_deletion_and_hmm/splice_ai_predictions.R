
setwd("/omics/groups/OE0433/internal/pablo/projects/sandbox/illumina/")

library(tidyverse)

all_exons_table <- read_delim(file = "data/209_GTEx_Exon_Regions.bed",
                              delim = "\t",
                              col_names = TRUE, guess_max = 99999)

col_types <- cols(
  CHROM = col_double(),
  POS = col_double(),
  ID = col_character(),
  REF = col_character(),
  ALT = col_character(),
  QUAL = col_character(),
  FILTER = col_character(),
  INFO = col_character()
)
all_exon_ids <- all_exons_table$ID

for (my_exon_id in all_exon_ids) {
  
  # my_exon_id <- "SE_10_+_90770357_90770510_90770572_90771756_FAS"
  
  print(which(all_exon_ids == my_exon_id))
  
  # from annotation
  strand <- strsplit(x = my_exon_id,
                     split = "_")[[1]][3]
  
  exon_coordinate_start <- as.numeric(strsplit(x = my_exon_id,
                                               split = "_")[[1]][5])
  exon_coordinate_end <- as.numeric(strsplit(x = my_exon_id,
                                             split = "_")[[1]][6])
  exon_chromosome <- strsplit(x = my_exon_id,
                              split = "_")[[1]][2]
  
  
  # initialise empty vectors that I'll fill in with information
  vector_exon_ids <- c()
  vector_deletion_ids <- c()
  vector_predictions <- c()
  
  
  if (strand == "+") {
    my_tabix_command <- paste("tabix data/spliceai_scores.raw.indel.hg19.vcf.gz ",
                              exon_chromosome,
                              ":",
                              exon_coordinate_start-1,
                              "-",
                              exon_coordinate_end,
                              " > data/tmp/",
                              gsub(pattern = "/", replacement = "-", x = my_exon_id),
                              ".tsv",
                              sep = "",
                              collapse = "")
    
    # run tabix command
    system(command = my_tabix_command)
    
    
    splice_ai_predictions <- read_delim(file = paste("data/tmp/",
                                                     gsub(pattern = "/", replacement = "-", x = my_exon_id),
                                                     ".tsv",
                                                     sep = "",
                                                     collapse = ""),
                                        delim = "\t",
                                        col_names = c("CHROM",
                                                      "POS",
                                                      "ID",
                                                      "REF",
                                                      "ALT",
                                                      "QUAL",
                                                      "FILTER",
                                                      "INFO"),
                                        col_types = col_types)
    
    if (nrow(splice_ai_predictions) == 0) {
      next
    }
    
    splice_ai_predictions <- splice_ai_predictions %>%
      filter(nchar(REF) > nchar(ALT)) # filter for deletions only
    
    exon_length <- exon_coordinate_end - exon_coordinate_start + 1
    
    
    # from VCF
    exon_position_minus_1 <- exon_coordinate_start - 1
    
    for (i in 1:exon_length){
      # i <- 10
      distance_to_exon_end <- exon_length - i
      
      this_coordinate <- exon_position_minus_1 + i - 1
      
      idx_this_position <- which(splice_ai_predictions$POS == this_coordinate)
      
      for (j in idx_this_position) {
        # j <- 48
        deletion_length <- nchar(splice_ai_predictions$REF[j]) - nchar(splice_ai_predictions$ALT[j])
        
        # skip this line if the deletion length is longer than the remainder of the exon
        if (deletion_length > distance_to_exon_end + 1) {
          next
        }
        
        
        this_deletion <- paste("Del_k",
                               deletion_length,
                               "_",
                               i,
                               "-",
                               i+deletion_length-1,
                               sep = "",
                               collapse = "")
        
        acceptor_gain_score <- as.numeric(strsplit(x = splice_ai_predictions$INFO[j],
                                                   split = "|",
                                                   fixed = TRUE)[[1]][3])
        acceptor_loss_score <- -1*as.numeric(strsplit(x = splice_ai_predictions$INFO[j],
                                                      split = "|",
                                                      fixed = TRUE)[[1]][4])
        donor_gain_score <- as.numeric(strsplit(x = splice_ai_predictions$INFO[j],
                                                split = "|",
                                                fixed = TRUE)[[1]][5])
        donor_loss_score <- -1*as.numeric(strsplit(x = splice_ai_predictions$INFO[j],
                                                   split = "|",
                                                   fixed = TRUE)[[1]][6])
        
        prediction <- c(acceptor_gain_score, acceptor_loss_score, donor_gain_score, donor_loss_score)[which(abs(c(acceptor_gain_score, acceptor_loss_score, donor_gain_score, donor_loss_score)) == max(abs(c(acceptor_gain_score, acceptor_loss_score, donor_gain_score, donor_loss_score))))]
        # print(this_deletion)
        # print(prediction)
        prediction <- prediction[1]
        
        vector_exon_ids <- c(vector_exon_ids,
                             my_exon_id)
        vector_deletion_ids <- c(vector_deletion_ids,
                                 this_deletion)
        vector_predictions <- c(vector_predictions,
                                prediction)
        
      }
    }
  } else if (strand == "-") {
    my_tabix_command <- paste("tabix data/spliceai_scores.raw.indel.hg19.vcf.gz ",
                              exon_chromosome,
                              ":",
                              exon_coordinate_end-1,
                              "-",
                              exon_coordinate_start,
                              " > data/tmp/",
                              gsub(pattern = "/", replacement = "-", x = my_exon_id),
                              ".tsv",
                              sep = "",
                              collapse = "")
    
    # run tabix command
    system(command = my_tabix_command)
    
    
    splice_ai_predictions <- read_delim(file = paste("data/tmp/",
                                                     gsub(pattern = "/", replacement = "-", x = my_exon_id),
                                                     ".tsv",
                                                     sep = "",
                                                     collapse = ""),
                                        delim = "\t",
                                        col_names = c("CHROM",
                                                      "POS",
                                                      "ID",
                                                      "REF",
                                                      "ALT",
                                                      "QUAL",
                                                      "FILTER",
                                                      "INFO"),
                                        col_types = col_types) 
    
    if (nrow(splice_ai_predictions) == 0) {
      next
    }
    
    splice_ai_predictions <- splice_ai_predictions %>%
      filter(nchar(REF) > nchar(ALT)) # filter for deletions only
    
    exon_length <- exon_coordinate_start - exon_coordinate_end + 1
    
    
    # from VCF
    exon_position_minus_1 <- exon_coordinate_start + 1
    
    vector_exon_ids <- c()
    vector_deletion_ids <- c()
    vector_predictions <- c()
    
    
    for (i in exon_length:1){
      # i <- 157
      distance_to_exon_end <- exon_length - i
      
      this_coordinate <- exon_position_minus_1 - i - 1
      
      idx_this_position <- which(splice_ai_predictions$POS == this_coordinate)
      
      for (j in idx_this_position) {
        # j <- 22
        deletion_length <- nchar(splice_ai_predictions$REF[j]) - nchar(splice_ai_predictions$ALT[j])
        
        # skip this line if the deletion length is longer than the remainder of the exon
        if (deletion_length > distance_to_exon_end + 1) {
          next
        }
        
        
        this_deletion <- paste("Del_k",
                               deletion_length,
                               "_",
                               i,
                               "-",
                               i+deletion_length-1,
                               sep = "",
                               collapse = "")
        
        acceptor_gain_score <- as.numeric(strsplit(x = splice_ai_predictions$INFO[j],
                                                   split = "|",
                                                   fixed = TRUE)[[1]][3])
        acceptor_loss_score <- -1*as.numeric(strsplit(x = splice_ai_predictions$INFO[j],
                                                      split = "|",
                                                      fixed = TRUE)[[1]][4])
        donor_gain_score <- as.numeric(strsplit(x = splice_ai_predictions$INFO[j],
                                                split = "|",
                                                fixed = TRUE)[[1]][5])
        donor_loss_score <- -1*as.numeric(strsplit(x = splice_ai_predictions$INFO[j],
                                                   split = "|",
                                                   fixed = TRUE)[[1]][6])
        
        prediction <- c(acceptor_gain_score, acceptor_loss_score, donor_gain_score, donor_loss_score)[which(abs(c(acceptor_gain_score, acceptor_loss_score, donor_gain_score, donor_loss_score)) == max(abs(c(acceptor_gain_score, acceptor_loss_score, donor_gain_score, donor_loss_score))))]
        # print(this_deletion)
        # print(prediction)
        prediction <- prediction[1]
        
        vector_exon_ids <- c(vector_exon_ids,
                             my_exon_id)
        vector_deletion_ids <- c(vector_deletion_ids,
                                 this_deletion)
        vector_predictions <- c(vector_predictions,
                                prediction)
        
      }
    }
  }
  
  
  write_delim(x = tibble(exon = vector_exon_ids,
                         exon_length = exon_length,
                         mutation = vector_deletion_ids,
                         prediction = vector_predictions),
              file = paste("data/exon_predictions/",
                           gsub(pattern = "/", replacement = "-", x = my_exon_id),
                           ".tsv",
                           sep = "",
                           collapse = ""),
              delim = "\t",
              col_names = TRUE, )
  
} # finish looping through all exons

# plot(vector_predictions[grep(pattern = "k4", x = vector_deletion_ids)], type = "l")
# abline(h=0)






















