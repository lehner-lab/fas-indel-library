
library(tidyverse)
library(scales)
load("data/indels_library.RData")


wt_sequence <- indels_library$sequence[which(indels_library$id == "WT")]
idx_singles <- which(nchar(indels_library$sequence)==63 & indels_library$id != "WT")
density_object <- density(x = indels_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]

indels_library$es <- indels_library$es - library_bias
indels_library$es[which(indels_library$id == "WT")] <- 0

indels_library$psi <- (49.1/1) * exp(indels_library$es)

indels_library$psi[which(indels_library$psi>100)] <-  100




indels_library_position_40 <- indels_library %>%
  filter(grepl(pattern = "_40_", x = id))

idx_ag <- c(grep(pattern = "_40_[TC]*AG",
                 x = indels_library_position_40$id),
            grep(pattern = "_40_[TCAG]*A$",
                 x = indels_library_position_40$id),
            grep(pattern = "_40_[TCAG]*A;",
                 x = indels_library_position_40$id))

indels_library_position_40$AG <- "no"
indels_library_position_40$AG[idx_ag] <- "yes"

table(indels_library_position_40$AG )

ggplot(data = indels_library_position_40,
       mapping = aes(x = factor(AG, levels = c("yes", "no")),
                     y = psi,
                     group = factor(AG, levels = c("yes", "no")))) +
  geom_violin(fill = "gray90", colour = NA, scale = "width") +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1.5) +
  ylab("percent spliced in") +
  xlab("creates an AG?") +
  ggpubr::stat_compare_means(comparisons = list(c("yes", "no")),
                             tip.length = 0,
                             method = "wilcox.test",
                             step.increase = 0.01)


ggsave(filename = "116_insert_AG_position_40.pdf",
       height = 5,
       width = 5,
       useDingbats = F)


write_delim(x = indels_library_position_40,
            delim = "\t",
            file = "116_figure_03o.tsv")
