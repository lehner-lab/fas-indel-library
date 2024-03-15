
setwd("/omics/groups/OE0433/internal/pablo/projects/sandbox/illumina/")

library(tidyverse)



load("hmm_tibble_exome_overview.RData")


ggplot(data = tibble_exome_overview %>%
         filter(length <= 200 & length >= 50),
       mapping = aes(x = length,
                     y =  sum_E)) +
  geom_point() +
  geom_smooth(method = "loess", span = 1 ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.75)




ggplot(data = tibble_exome_overview %>%
         filter(length <= 200 & length >= 50),
       mapping = aes(x = length,
                     y =  sum_S+sum_E)) +
  geom_point() +
  geom_smooth(method = "loess", span = 5 ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.75)


y50 <- mean((tibble_exome_overview$sum_E + tibble_exome_overview$sum_S)[tibble_exome_overview$length == 50])
y51 <-  mean((tibble_exome_overview$sum_E + tibble_exome_overview$sum_S)[tibble_exome_overview$length == 51])
y200 <- mean((tibble_exome_overview$sum_E + tibble_exome_overview$sum_S)[tibble_exome_overview$length == 200])

# y = mx + c
# 
# when x = 50, y = y50

# 
# y50 = 50m + c

# 
# c = y50 - 50m


m_ = y51/51

c_ = y51 - 51*m_






ggplot(data = tibble_exome_overview %>%
         filter(length <= 200 & length >= 50),
       mapping = aes(x = length,
                     y =  sum_S+sum_E)) +
  geom_point() +
  geom_smooth() +
  geom_abline(intercept = c_,
              slope = m_,
              colour = "red",
              lwd = 2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 0.75) +
  xlab("exon length") +
  ylab("sum of enhancer + silencer 4mers")
ggsave(filename = "plot_sum_of_E_S_vs_exon_length.pdf",
       useDingbats = F,
       width = 7,
       height = 7)