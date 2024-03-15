
library(tidyverse)


load("/Users/p028n/CRG/Indel/data/PB_FAS_Indel_SE_fitness_replicates.RData")


correlation_table <- all_variants[,24:32]
colnames(correlation_table) <- paste("rep ", 1:9, sep = "")
# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y,use = "complete.obs"), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}

pdf(file = "142_correlation_plot.pdf", width = 8, height = 8, useDingbats = F)
# Create the plots
pairs(all_variants[,24:32], 
      lower.panel = panel.cor,
      upper.panel = upper.panel)
dev.off()


