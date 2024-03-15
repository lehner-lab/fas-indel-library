setwd("/omics/groups/OE0433/internal/pablo/projects/sandbox/bp_ppt")


# load table
Bed.Table <- read.table(file = "/omics/groups/OE0433/internal/pablo/projects/sandbox/bp_ppt/bed_file_with_mRNA_sequences.bed",
                        header = T,
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




list_bp_scores <- vector(mode = "list",
                         length = 0)



for (each_row in 1:nrow(Bed.Table)) {
  
  print(each_row)
  
  exon_name <- Bed.Table$id[each_row]
  exon_length <- Bed.Table$length[each_row]
  exon_sequence <- Bed.Table$sequence[each_row]
  
  
  
  #### MAKE FASTA FILE ####
  
  my_fasta_file <- data.frame(only_col = c())
  
  if (nchar(exon_sequence) >= 20) {
    
    sequence_length <- nchar(exon_sequence)
    
    for (i in 20:sequence_length) {
      
      sequence_start <- 1
      sequence_end <- i
      
      if (sequence_end > 500) {
        sequence_start <- sequence_end - 500 + 1
      }
      
      my_fasta_file <- rbind(my_fasta_file,
                             data.frame(only_col = c(paste(">",
                                                           exon_name,
                                                           ":",
                                                           i,
                                                           sep = "",
                                                           collapse = ""),
                                                     substr(x = exon_sequence,
                                                            start = sequence_start,
                                                            stop = sequence_end))))
    }
    
  }
  
  
  fasta_filename <- paste("/omics/groups/OE0433/internal/pablo/projects/sandbox/bp_ppt/fasta/",
                          each_row,
                          ".fa",
                          sep = "",
                          collapse = "")
  
  write.table(x = my_fasta_file,
              file = fasta_filename,
              quote = F,
              col.names = F,
              row.names = F,
              sep = "\t")
  
  
  
  
  
  
  
  
  
  
  #### RUN BP FINDER ####
  
  bpfinder_command <- paste("python2 ",
                            "/omics/groups/OE0433/internal/pablo/projects/sandbox/bp_ppt/regulatorygenomicsupf-svm-bpfinder-727e2d8ea4ad/svm_bpfinder.py ",
                            "-i /omics/groups/OE0433/internal/pablo/projects/sandbox/bp_ppt/fasta/", each_row, ".fa ",
                            "-s Hsap ",
                            "-l 1000 ",
                            "-d 15 ",
                            "> ",
                            "/omics/groups/OE0433/internal/pablo/projects/sandbox/bp_ppt/bp_finder_output/", each_row, ".output",
                            sep = "",
                            collapse = ""
                            )
  
  system(bpfinder_command)
  
  
  
  
  
  
  
  
  #### GET MAX SCORE FROM BP FINDER
  
  bpfinder_output <- read.table(file = paste("/omics/groups/OE0433/internal/pablo/projects/sandbox/bp_ppt/bp_finder_output/",
                                             each_row,
                                             ".output",
                                             sep = "",
                                             collapse = ""),
                                header = T,
                                sep = "\t")
  
  idx_max <- which(bpfinder_output$svm_scr == max(bpfinder_output$svm_scr))[1]
  
  if (is.na(idx_max)) {
    max_bp_score <- as.data.frame(matrix(-999, nrow = 1, ncol = 10))
  } else {
    max_bp_score <- bpfinder_output[idx_max,]
  }
  
  
  
  
  list_bp_scores[[length(list_bp_scores)+1]] <- max_bp_score
  
  
}




bp_scores_df <- Reduce(rbind, list_bp_scores)
bp_scores_df <- data.table::rbindlist(list_bp_scores)

save(bp_scores_df, file = "bp_scores.RData")


write.table(bp_scores_df,
            file = "bp_scores.txt",
            sep = "\t",
            quote = F,
            col.names = T,
            row.names = F)
save(bp_scores_df, Bed.Table, file = "bp_scores.RData")

mean(bp_scores_df$svm_scr > 1.1901840)*100


library(tidyverse)







# 
# hist(bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)])
# hist(bp_scores_df$bp_scr[which(bp_scores_df$svm_scr>-900)])
# hist(bp_scores_df$ppt_scr[which(bp_scores_df$svm_scr>-900)])
# 
# View(bp_scores_df[which(Bed.Table$gene=="FAS"),])
# plot(log2(Bed.Table$length[which(bp_scores_df$svm_scr>-900)]),
#      bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)],
#      pch = ".")
# plot(log2(Bed.Table$length[which(bp_scores_df$svm_scr>-900)]),
#      bp_scores_df$bp_scr[which(bp_scores_df$svm_scr>-900)],
#      pch = ".")
# plot(log2(Bed.Table$length[which(bp_scores_df$svm_scr>-900)]),
#      log2(bp_scores_df$ppt_scr[which(bp_scores_df$svm_scr>-900)]))
# plot(log2(Bed.Table$length[which(bp_scores_df$svm_scr>-900)]),
#      log2(bp_scores_df$ppt_len[which(bp_scores_df$svm_scr>-900)]))
# 
# plot(bp_scores_df$ppt_scr[which(bp_scores_df$svm_scr>-900)],
#      log2(bp_scores_df$ppt_len[which(bp_scores_df$svm_scr>-900)]))
# 
# plot(bp_scores_df$bp_scr[which(bp_scores_df$svm_scr>-900)],
#      log2(bp_scores_df$ppt_scr[which(bp_scores_df$svm_scr>-900)]))
# 
# 
# boxplot(bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] ~ findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                                                               vec = seq(0,500,25), all.inside = T))
# hist(log2(Bed.Table$length[which(bp_scores_df$svm_scr>-900)]))
# 
# 
# boxplot(bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] ~ Bed.Table$strand[which(bp_scores_df$svm_scr>-900)])
# 
# 
# mean(bp_scores_df$svm_scr > 1)*100

ggplot(data = bp_scores_df[which(bp_scores_df$svm_scr>-900),],
       mapping = aes(x = svm_scr)) +
  geom_density() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab("SVM score") +
  geom_vline(xintercept = 1.1901840,
             lty = 2, 
             colour = "red")
# 
# mean(bp_scores_df$svm_scr[Bed.Table$length < 70] > 1.1901840)*100
# 
# mean(bp_scores_df$bp_scr > 1)*100
# mean(bp_scores_df$bp_scr > 2.1400875)*100
# 
# 
# mean(bp_scores_df$ppt_scr > 44)*100
# mean(bp_scores_df$ppt_scr > 44 & bp_scores_df$bp_scr > 2.1400875)*100
# 


load("psi_values/001_PSI_Values.RData")
Many.NAs <- which(apply(X = Merged.PSI[,2:ncol(Merged.PSI)],
                        MARGIN = 1,
                        FUN = function(x){
                          mean(is.na(x)) > 0.8
                        }))
Merged.PSI <- Merged.PSI[-Many.NAs,]

Exon.Lengths <- sapply(X = as.character(Merged.PSI$rn),
                       FUN = function(x){
                         coordinates <- as.numeric(strsplit(x,"_")[[1]][5:6])
                         exon.length <- abs(coordinates[2] - coordinates[1]) + 1
                         exon.length
                       })
Merged.PSI <- Merged.PSI[-which(Exon.Lengths<2),]


count <- 1
mean_psi <- sapply(X = as.character(Bed.Table$id),
                   FUN = function(x){
                     print(count)
                     count <<- count + 1
                     idx <- which(as.character(Merged.PSI$rn) == x)
                     mean(as.numeric(Merged.PSI[idx,2:ncol(Merged.PSI)]), na.rm = T)
                   })

# mean(as.numeric(Merged.PSI[13068,2:ncol(Merged.PSI)]), na.rm = T)
# mean(Merged.PSI[13068,2:10], na.rm = T)
# 
# 
# plot(bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)],
#      mean_psi[which(bp_scores_df$svm_scr>-900)])




# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] > 1.1901840,
#                      psi = 100* mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                            vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(x = as.factor(length),
#                      y = psi,
#                      fill = svm_score)) +
#   geom_boxplot(notch = T) +
#   scale_fill_manual(values = c("firebrick1", "orange"),
#                     name = "svm score > 1") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_x_discrete(breaks=c("1","2","3","4","5"),
#                    labels=c("[0 – 100)", "[100 – 200)", "[200 – 300)", "[300 – 400)", "[400 – 500]")) +
#   xlab("exon length") +
#   ylab("PSI")


ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] > 1.1901840,
                     psi = mean_psi[which(bp_scores_df$svm_scr>-900)] * 100,
                     length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
                                           vec = seq(0,500,100), all.inside = T)),
       mapping = aes(x = as.factor(length),
                     y = psi,
                     fill = svm_score)) +
  geom_violin(scale = "width",
              position = position_dodge(width=0.9)) +
  stat_summary(mapping = aes(x = as.factor(length),
                             y = psi,
                             fill = svm_score),
               geom = "point",
               position = position_dodge(width=0.9),
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +
  scale_fill_manual(values = c("#ECECEC", "#F9C45F"),
                    name = "svm score > FAS exon 6") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.85) +
  scale_x_discrete(breaks=c("1","2","3","4","5"),
                   labels=c("[0 – 100)", "[100 – 200)", "[200 – 300)", "[300 – 400)", "[400 – 500]")) +
  xlab("exon length") +
  ylab("PSI")


ggsave(filename = "psi_vs_SVM_vs_exon_lengthplot.pdf", width = 8, height = 8, useDingbats = F)




# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] > 1,
#                      psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      length = Bed.Table$length[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(x = length,
#                      y = psi,
#                      fill = as.factor(svm_score),
#                      colour = as.factor(svm_score),
#                      group = as.factor(svm_score))) +
#   geom_smooth(method = "loess", span = 1) +
#   xlim(c(0,500)) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_fill_manual(values = c("firebrick1", "orange"),
#                     name = "svm score > 1") +
#   scale_colour_manual(values = c("firebrick1", "orange"),
#                     name = "svm score > 1") +
#   xlab("exon length") +
#   ylab("PSI")

# 
# ggplot(data = tibble(svm_score = findInterval(x = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)],
#                                               vec = seq(0,2,0.5), all.inside = T),
#                      psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      length = Bed.Table$length[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(x = length,
#                      y = psi,
#                      fill = as.factor(svm_score),
#                      group = as.factor(svm_score))) +
#   geom_smooth(method = "loess", span = 1) +
#   coord_cartesian(xlim = c(0,500))
# 
# 
# 
# 
# 
# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)],
#                      length = log2(Bed.Table$length[which(bp_scores_df$svm_scr>-900)])),
#        mapping = aes(x = length,
#                      y = svm_score)) +
#   geom_point(alpha = 0.25) +
#   geom_smooth(fill = "red",
#               colour = "red") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("log2 exon length") +
#   ylab("SVM score")
# 
# 
# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)],
#                      psi = mean_psi[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(x = svm_score,
#                      y = psi)) +
#   geom_point(alpha = 0.25) +
#   geom_smooth(fill = "red",
#               colour = "red") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("SVM score") +
#   ylab("PSI")
# 

# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)],
#                      psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                            vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(x = svm_score,
#                      y = psi,
#                      group = as.factor(length),
#                      colour = as.factor(length),
#                      fill = as.factor(length))) +
#   geom_point(alpha = 0.25,
#              colour = "black") +
#   geom_smooth() +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("SVM score") +
#   ylab("PSI")





# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] > 2,
#                      psi = 100* mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                            vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(x = as.factor(length),
#                      y = psi,
#                      fill = svm_score)) +
#   geom_boxplot(notch = T) +
#   scale_fill_manual(values = c("firebrick1", "orange"),
#                     name = "svm score > 2") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_x_discrete(breaks=c("1","2","3","4","5"),
#                    labels=c("[0 – 100)", "[100 – 200)", "[200 – 300)", "[300 – 400)", "[400 – 500]")) +
#   xlab("exon length") +
#   ylab("PSI")
# 
# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] > 2,
#                      psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                            vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(x = as.factor(length),
#                      y = psi,
#                      fill = svm_score)) +
#   geom_violin(scale = "width") +
#   scale_fill_manual(values = c("firebrick1", "orange"),
#                     name = "svm score > 2") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_x_discrete(breaks=c("1","2","3","4","5"),
#                    labels=c("[0 – 100)", "[100 – 200)", "[200 – 300)", "[300 – 400)", "[400 – 500]")) +
#   xlab("exon length") +
#   ylab("PSI")






# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] > 2,
#                      psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      length = Bed.Table$length[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(x = length,
#                      y = psi,
#                      fill = as.factor(svm_score),
#                      colour = as.factor(svm_score),
#                      group = as.factor(svm_score))) +
#   geom_smooth(method = "loess", span = 1) +
#   xlim(c(0,500)) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_fill_manual(values = c("firebrick1", "orange"),
#                     name = "svm score > 2") +
#   scale_colour_manual(values = c("firebrick1", "orange"),
#                       name = "svm score > 2") +
#   xlab("exon length") +
#   ylab("PSI")
# 





# https://rdrr.io/cran/Peptides/man/hydrophobicity.html
# https://rdrr.io/cran/seqinr/man/translate.html


count <- 0
hydrophobicity_scores <- sapply(X = Bed.Table$sequence,
                                FUN = function(x){
                                  count <<- count+1
                                  print(count)
                                  
                                  final_hydrophobicity_score <- NA
                                  
                                  if (nchar(x) >= 5) {
                                    frame_0 <- seqinr::translate(seq = strsplit(x = x,
                                                                                split = "")[[1]],
                                                                 frame = 0)
                                    
                                    
                                    frame_1 <- seqinr::translate(seq = strsplit(x = x,
                                                                                split = "")[[1]],
                                                                 frame = 1)
                                    
                                    
                                    frame_2 <- seqinr::translate(seq = strsplit(x = x,
                                                                                split = "")[[1]],
                                                                 frame = 2)
                                    
                                    
                                    
                                    hydrophobicity_0 <- NA
                                    hydrophobicity_1 <- NA
                                    hydrophobicity_2 <- NA
                                    
                                    
                                    
                                    if (! "*" %in% frame_0) {
                                      hydrophobicity_0 <- Peptides::hydrophobicity(seq = paste(frame_0,
                                                                                               sep = "",
                                                                                               collapse = ""),
                                                                                   scale = "KyteDoolittle")
                                    }
                                    
                                    if (! "*" %in% frame_1) {
                                      hydrophobicity_1 <- Peptides::hydrophobicity(seq = paste(frame_1,
                                                                                               sep = "",
                                                                                               collapse = ""),
                                                                                   scale = "KyteDoolittle")
                                    }
                                    
                                    if (! "*" %in% frame_2) {
                                      hydrophobicity_2 <- Peptides::hydrophobicity(seq = paste(frame_2,
                                                                                               sep = "",
                                                                                               collapse = ""),
                                                                                   scale = "KyteDoolittle")
                                    }
                                    
                                    
                                    
                                    
                                    if (! all(is.na(c(hydrophobicity_0, hydrophobicity_1, hydrophobicity_2)))) {
                                      final_hydrophobicity_score <- max(c(hydrophobicity_0, hydrophobicity_1, hydrophobicity_2),
                                                                        na.rm = T)
                                    }
                                  }
                                  
                                  final_hydrophobicity_score
                                  
                                })

# 
# 
# min(c(NA,NA), na.rm = T)
# Peptides::hydrophobicity("WWWW", scale = "Argos")
# Peptides::hydrophobicity("AAAA", scale = "Argos")
# Peptides::hydrophobicity("QQH", scale = "Argos")
# 
# Peptides::hydrophobicity("WWWW", scale = "Juretic")
# Peptides::hydrophobicity("AAAAW", scale = "Juretic")
# Peptides::hydrophobicity("QQH", scale = "Juretic")
# 
# 
# Peptides::hydrophobicity("WWWW", scale = "KyteDoolittle")
# Peptides::hydrophobicity("AAAAW", scale = "KyteDoolittle")
# Peptides::hydrophobicity("QQH", scale = "KyteDoolittle")
# 








# 
# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(x = svm_score,
#                      y = hydropathy)) +
#   geom_point(alpha = 0.25) +
#   geom_smooth(fill = "red",
#               colour = "red") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("SVM score") +
#   ylab("Hydropathy")


# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                            vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(x = svm_score,
#                      y = hydropathy,
#                      group = as.factor(length),
#                      colour = as.factor(length),
#                      fill = as.factor(length))) +
#   geom_point(alpha = 0.25,
#              colour = "black") +
#   geom_smooth() +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("SVM score") +
#   ylab("Hydropthy") 






# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                            vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(y = psi,
#                      x = hydropathy,
#                      group = as.factor(length),
#                      colour = as.factor(length),
#                      fill = as.factor(length))) +
#   geom_point(alpha = 0.25,
#              colour = "black") +
#   geom_smooth(method = "loess", span = 2) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("Hydropthy") +
#   ylab("PSI") 










# 
# ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] > 1,
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                            vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(x = as.factor(length),
#                      y = hydropathy,
#                      fill = svm_score)) +
#   geom_boxplot(notch = T) +
#   scale_fill_manual(values = c("firebrick1", "orange"),
#                     name = "svm score > 1") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   scale_x_discrete(breaks=c("1","2","3","4","5"),
#                    labels=c("[0 – 100)", "[100 – 200)", "[200 – 300)", "[300 – 400)", "[400 – 500]")) +
#   xlab("exon length") +
#   ylab("Hydropathy")


ggplot(data = tibble(svm_score = bp_scores_df$svm_scr[which(bp_scores_df$svm_scr>-900)] > 1.1901840,
                     hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
                     length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
                                           vec = seq(0,500,100), all.inside = T)),
       mapping = aes(x = as.factor(length),
                     y = hydropathy,
                     fill = svm_score)) +
  geom_hline(yintercept = 0,
             lty = 3) +
  geom_violin(scale = "width",
              position = position_dodge(width=0.9)) +
  stat_summary(mapping = aes(x = as.factor(length),
                             y = hydropathy,
                             fill = svm_score),
               geom = "point",
               position = position_dodge(width=0.9),
               # shape = "-",
               size = 3,
               colour = "black",
               # fun.min = function(z) { quantile(z,0.25) },
               # fun.max = function(z) { quantile(z,0.75) },
               fun = median) +
  scale_fill_manual(values = c("#ECECEC", "#F9C45F"),
                    name = "svm score > FAS exon 6") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.85) +
  scale_x_discrete(breaks=c("1","2","3","4","5"),
                   labels=c("[0 – 100)", "[100 – 200)", "[200 – 300)", "[300 – 400)", "[400 – 500]")) +
  xlab("exon length") +
  ylab("Hydropathy")


ggsave(filename = "hydropathy_plot.pdf", width = 8, height = 8, useDingbats = F)





# 
# 
# 
# 
# 
# distance_from_3ss_to_end <- apply(X = tibble(id = bp_scores_df$seq_id,
#                                              dist_bp_to_3ss = bp_scores_df$ss_dist,
#                                              exon_length = Bed.Table$length),
#                                   MARGIN = 1,
#                                   FUN = function(x){
#                                     sequence_id <- as.character(x[1])
#                                     distance_bp_3ss <- as.numeric(x[2])
#                                     exon_length <- as.numeric(x[3])
#                                     
#                                     position_3ss_in_exon <- as.numeric(strsplit(sequence_id, ":")[[1]][2])
#                                     
#                                     distance_from_3ss_to_end <- exon_length - position_3ss_in_exon
#                                   })
# 
# distance_from_3ss_to_beginning <- apply(X = tibble(id = bp_scores_df$seq_id,
#                                                    dist_bp_to_3ss = bp_scores_df$ss_dist,
#                                                    exon_length = Bed.Table$length),
#                                         MARGIN = 1,
#                                         FUN = function(x){
#                                           sequence_id <- as.character(x[1])
#                                           distance_bp_3ss <- as.numeric(x[2])
#                                           exon_length <- as.numeric(x[3])
#                                           
#                                           position_3ss_in_exon <- as.numeric(strsplit(sequence_id, ":")[[1]][2])
#                                           
#                                           position_3ss_in_exon
#                                         })
# 
# distance_from_bp_to_end <- apply(X = tibble(id = bp_scores_df$seq_id,
#                                             dist_bp_to_3ss = bp_scores_df$ss_dist,
#                                             exon_length = Bed.Table$length),
#                                  MARGIN = 1,
#                                  FUN = function(x){
#                                    sequence_id <- as.character(x[1])
#                                    distance_bp_3ss <- as.numeric(x[2])
#                                    exon_length <- as.numeric(x[3])
#                                    
#                                    position_3ss_in_exon <- as.numeric(strsplit(sequence_id, ":")[[1]][2])
#                                    
#                                    exon_length - (position_3ss_in_exon - distance_bp_3ss)
#                                  })
# 
# distance_from_bp_to_beginning <- apply(X = tibble(id = bp_scores_df$seq_id,
#                                                   dist_bp_to_3ss = bp_scores_df$ss_dist,
#                                                   exon_length = Bed.Table$length),
#                                        MARGIN = 1,
#                                        FUN = function(x){
#                                          sequence_id <- as.character(x[1])
#                                          distance_bp_3ss <- as.numeric(x[2])
#                                          exon_length <- as.numeric(x[3])
#                                          
#                                          position_3ss_in_exon <- as.numeric(strsplit(sequence_id, ":")[[1]][2])
#                                          
#                                          position_3ss_in_exon - distance_bp_3ss
#                                        })
# 
# 
# # View(Bed.Table[which(Bed.Table$gene == "FAS"),])
# # View(bp_scores_df[which(Bed.Table$gene == "FAS"),])
# # 
# # 
# # strsplit("SE_10_+_90770357_90770510_90770572_90771756_FAS:47", split = ":")
# # 
# # GATCCAGATC TAACTTGGGG TGGCTTTGTC TTCTTCTTTT GCCAATTCCA CTAATTGTTT GGG
# # 
# # 
# # GATCCAGATC TA ACTTGGGG TGGCTTTGTC TTCTTCTTTT GCCAATT CCA CTAATTGTTT GGG
# # 
# # #distance from new 3'ss to the end
# # 63-47
# # # distance from new 3'ss to beginning
# # 47
# # 
# # # distance from bp to the end
# # 
# # 63 - (47-35)
# # 
# # 
# # # distance from bp to the beginning
# # 47-35
# 
# 
# hist(log2(distance_from_3ss_to_end))
# 
# 
# 
# 
# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = log2(distance_from_bp_to_beginning)[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(x = length#,
#                      # group = as.factor(length),
#                      # colour = as.factor(length),
#                      # fill = as.factor(length)
#        )) +
#   geom_density() +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("log2 distance from BP to start of exon") 
# 
# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = log2(distance_from_3ss_to_end)[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(x = length#,
#                      # group = as.factor(length),
#                      # colour = as.factor(length),
#                      # fill = as.factor(length)
#        )) +
#   geom_density() +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("log2 distance from py tract to end of exon") 
# 
# 
# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = log2(distance_from_bp_to_beginning)[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(y = psi,
#                      x = length#,
#                      # group = as.factor(length),
#                      # colour = as.factor(length),
#                      # fill = as.factor(length)
#        )) +
#   geom_point(alpha = 0.25,
#              colour = "black") +
#   geom_smooth(fill = "red", colour = "red") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("log2 distance from BP to start of exon") +
#   ylab("PSI")
# 
# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = log2(distance_from_bp_to_beginning)[which(bp_scores_df$svm_scr>-900)],
#                      exon_length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                                 vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(y = psi,
#                      x = length,
#                      group = as.factor(exon_length),
#                      colour = as.factor(exon_length),
#                      fill = as.factor(exon_length)
#        )) +
#   geom_point(alpha = 0.25,
#              colour = "black") +
#   geom_smooth() +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("log2 distance from BP to start of exon") +
#   ylab("PSI")
# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = log2(distance_from_3ss_to_end)[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(y = psi,
#                      x = length#,
#                      # group = as.factor(length),
#                      # colour = as.factor(length),
#                      # fill = as.factor(length)
#        )) +
#   geom_point(alpha = 0.25,
#              colour = "black") +
#   geom_smooth(fill = "red", colour = "red") +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("log2 distance from py tract to end of exon") +
#   ylab("PSI")
# 
# 
# 
# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = log2(distance_from_3ss_to_end)[which(bp_scores_df$svm_scr>-900)],
#                      exon_length = findInterval(x = Bed.Table$length[which(bp_scores_df$svm_scr>-900)],
#                                                 vec = seq(0,500,100), all.inside = T)),
#        mapping = aes(y = psi,
#                      x = length,
#                      group = as.factor(exon_length),
#                      colour = as.factor(exon_length),
#                      fill = as.factor(exon_length)
#        )) +
#   geom_point(alpha = 0.25,
#              colour = "black") +
#   geom_smooth() +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("log2 distance from py tract to end of exon") +
#   ylab("PSI")
# 
# 
# 
# # ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
# #                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
# #                      length = (distance_from_3ss_to_end < 16)[which(bp_scores_df$svm_scr>-900)]),
# #        mapping = aes(y = psi,
# #                      x = length,
# #                      group = as.factor(length),
# #                      #colour = as.factor(length),
# #                      fill = as.factor(length)
# #        )) +
# #   geom_boxplot(notch = T)
# # 
# # 
# # hist(log2(distance_from_bp_to_beginning))
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = (findInterval(x = log2(distance_from_bp_to_beginning), vec = seq(2,10,2), all.inside = T))[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(y = psi,
#                      x = length,
#                      group = as.factor(length),
#                      #colour = as.factor(length),
#                      fill = as.factor(length)
#        )) +
#   geom_boxplot(notch = T) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("distance from BP to exon start") +
#   ylab("PSI")
# 
# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = (findInterval(x = log2(distance_from_bp_to_beginning), vec = seq(2,10,2), all.inside = T))[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(y = psi,
#                      x = length,
#                      group = as.factor(length),
#                      #colour = as.factor(length),
#                      fill = as.factor(length)
#        )) +
#   geom_violin(scale = "width")
# 
# 
# 
# hist(log2(distance_from_3ss_to_end))
# 
# 
# ggplot(data = tibble(psi = mean_psi[which(bp_scores_df$svm_scr>-900)],
#                      hydropathy = hydrophobicity_scores[which(bp_scores_df$svm_scr>-900)],
#                      length = (findInterval(x = log2(distance_from_3ss_to_end), vec = seq(2,10,2), all.inside = T))[which(bp_scores_df$svm_scr>-900)]),
#        mapping = aes(y = psi,
#                      x = length,
#                      group = as.factor(length),
#                      #colour = as.factor(length),
#                      fill = as.factor(length)
#        )) +
#   geom_boxplot(notch = T) +
#   theme_bw() +
#   theme(panel.grid = element_blank()) +
#   xlab("distance from py tract to exon end") +
#   ylab("PSI")
