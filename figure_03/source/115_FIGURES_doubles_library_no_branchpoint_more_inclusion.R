
library(tidyverse)
load("/Users/p028n/CRG/Doubles/data/doubles_library.RData")




wt_sequence <- doubles_library$sequence[which(doubles_library$id == "")]
idx_singles <- which( (! grepl(pattern = ";", x = doubles_library$id)) & doubles_library$id != "")
density_object <- density(x = doubles_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]

doubles_library$es <- doubles_library$es - library_bias
doubles_library$es[which(doubles_library$id == "")] <- 0

doubles_library$psi <- (49.1/1) * exp(doubles_library$es)

doubles_library$psi[which(doubles_library$psi>100)] <-  100





idx_no_bp <- intersect(grep(pattern = "12", x = doubles_library$id),
                       grep(pattern = "13", x = doubles_library$id))
idx_no_12 <- grep(pattern = "12", x = doubles_library$id)

idx_no_13 <- grep(pattern = "13", x = doubles_library$id)



doubles_library$branchpoint <- "CTAACT"
doubles_library$branchpoint[idx_no_12] <- "CTXACT"
doubles_library$branchpoint[idx_no_13] <- "CTAXCT"
doubles_library$branchpoint[idx_no_bp] <- "CTXXCT"



table(doubles_library$branchpoint)


ggplot(data = doubles_library,
       mapping = aes(x = factor(branchpoint,
                                levels = c("CTAACT",
                                           "CTAXCT",
                                           "CTXACT",
                                           "CTXXCT")),
                     y = psi,
                     group = factor(branchpoint,
                                    levels = c("CTAACT",
                                               "CTAXCT",
                                               "CTXACT",
                                               "CTXXCT")))) +
  geom_violin(fill = "gray90", colour = NA, scale = "width") +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  theme_bw() +
  ylab("percent spliced in") +
  xlab("") +
  ggpubr::stat_compare_means(comparisons = list(c("CTAACT",
                                                  "CTAXCT"),
                                                c("CTAXCT",
                                                  "CTXACT"),
                                                c("CTXACT",
                                                  "CTXXCT")),
                             tip.length = 0,
                             method = "wilcox.test",
                             step.increase = 0.01) +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1.5,
        axis.text.x = element_text(angle = 90)) +
  xlab("") +
  ylab("percent spliced in") +
  ggtitle("exonic branchpoint mutations") +
  coord_cartesian(ylim = c(0,110)) +
  scale_y_continuous(breaks = seq(0,100,25))



ggsave(filename = "115_branchpoint_mutations.pdf",
       height = 5,
       width = 5,
       useDingbats = F)



write_delim(x = doubles_library,
            delim = "\t",
            file = "115_figure_03n.tsv")






####################3333
plot_1 <- ggplot(data = doubles_library,
       mapping = aes(x = factor(branchpoint,
                                levels = c("CTAACT",
                                           "CTAXCT",
                                           "CTXACT",
                                           "CTXXCT")),
                     y = psi)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(comparisons = list(c("CTAACT",
                                                  "CTAXCT"),
                                                c("CTAXCT",
                                                  "CTXACT"),
                                                c("CTXACT",
                                                  "CTXXCT")),
                             tip.length = 0,
                             step.increase = 0.01,
                             method = "wilcox.test") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1.5,
        axis.text.x = element_text(angle = 90)) +
  xlab("") +
  ylab("percent spliced in") +
  ggtitle("exonic branchpoint mutations") +
  coord_cartesian(ylim = c(0,100))

plot_2 <- ggplot(data = doubles_library,
       mapping = aes(x = factor(branchpoint,
                                levels = c("CTAACT",
                                           "CTAXCT",
                                           "CTXACT",
                                           "CTXXCT")),
                     y = psi)) +
  geom_violin(scale = "width") +
  ggpubr::stat_compare_means(comparisons = list(c("CTAACT",
                                                  "CTAXCT"),
                                                c("CTAXCT",
                                                  "CTXACT"),
                                                c("CTXACT",
                                                  "CTXXCT")),
                             tip.length = 0,
                             step.increase = 0.01,
                             method = "wilcox.test") +
  geom_boxplot(width=0.15,
               outlier.colour = NA) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1.5,
        axis.text.x = element_text(angle = 90)) +
  xlab("") +
  ylab("percent spliced in") +
  ggtitle("exonic branchpoint mutations") +
  coord_cartesian(ylim = c(0,100))


ggsave(filename = "115_branchpoint_mutations.pdf",
       plot = ggpubr::ggarrange(plot_1,plot_2,ncol=2), width = 8.27, height = 5.83)









ggplot(data = doubles_library,
       mapping = aes(x = factor(branchpoint,
                                levels = c("CTAACT",
                                           "CTAXCT",
                                           "CTXACT",
                                           "CTXXCT")),
                     y = psi)) +
  ggdist::stat_halfeye(
    width = 1,
    # custom bandwidth
    #adjust = 0.81,
    # move geom to the right
    justification = -.25,
    # remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    # remove outliers
    outlier.colour = NA,
    width = 0.3
  ) +
  ggpubr::stat_compare_means(comparisons = list(c("CTAACT",
                                                  "CTAXCT"),
                                                c("CTAXCT",
                                                  "CTXACT"),
                                                c("CTXACT",
                                                  "CTXXCT")),
                             tip.length = 0,
                             step.increase = 0.01,
                             method = "wilcox.test") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1,
        axis.text.x = element_text(angle = 90)) +
  xlab("") +
  ylab("percent spliced in") +
  ggtitle("exonic branchpoint mutations") +
  coord_cartesian(ylim = c(0,100))


ggsave(filename = "115_branchpoint_mutations_v2.pdf", width = 8, height = 8)




