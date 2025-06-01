
library(tidyverse)
library(scales)
load("data/indels_library.RData")


wt_sequence <- indels_library$sequence[which(indels_library$id == "WT")]
idx_singles <- which(nchar(indels_library$sequence)==63 & indels_library$id != "WT")
density_object <- density(x = indels_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]

indels_library$es <- indels_library$es - library_bias

indels_library$psi <- (49.1/1) * exp(indels_library$es)

indels_library$psi[which(indels_library$psi>100)] <-  100








psi_vector <- c()
position_vector <- c()
k_vector <- c()

# sliding windows size 1-10
# each window deletion is defined by the one and only deletion that exists there
for (k in 1:10) {
  
  # each position
  for (p in 1:(63-k+1)) {
    
    this_id <- paste("Del_k",
                     k,
                     "_",
                     p,
                     "-",
                     p+k-1,
                     sep = "",
                     collapse = "")
    
    idx <- grep(pattern = this_id, x = indels_library$id)
    
    if (length(idx)>0) {
      psi_vector <- c(psi_vector,
                      indels_library$psi[idx])
      position_vector <- c(position_vector,
                           p + (k/2) - 0.5)
      k_vector <- c(k_vector,
                    k)
    }
    
    
  }
  
  
}





plot_tibble <- tibble(psi = psi_vector,
                      k = k_vector,
                      position = position_vector)


ggplot(data = plot_tibble,
       mapping = aes(x = position,
                     y = psi,
                     colour = as.factor(k))) +
  geom_hline(yintercept = 49.1,
             colour = "gray80",
             size = 2) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank())



colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 11, palette = "Spectral")[2:10]))
my_colours <- colour_function(11)


ggplot(data = plot_tibble %>% filter(k %in% 1:6),
       mapping = aes(x = position,
                     y = psi)) +
  geom_hline(yintercept = 49.1,
             colour = "gray90",
             size = 2) +
  geom_point(colour = "gray30",
             fill = "gray30",
             mapping = aes(shape = as.factor(k))) +
  geom_smooth(method = "loess", span = 0.15, colour = my_colours[3], fill = my_colours[3], alpha = 0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 20/63) +
  scale_x_continuous(limits = c(1,63),
                     expand = c(0.01, 0.01),
                     breaks = 1:63) +
  scale_shape_manual(values = c(21:25,3),
                     name = "deletion length")
  
ggsave(filename = "127_exon_modules_plot.pdf",
       width = 10,
       height = 10,
       useDingbats = F)



write_delim(x = plot_tibble %>% filter(k %in% 1:6),
            delim = "\t",
            file = "127_figure_05c_dms.tsv")

# module 1

pdf(file = "127_module_1.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 9 & plot_tibble$position >= 1] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()

pdf(file = "127_module_2.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 18 & plot_tibble$position >= 10] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()

pdf(file = "127_module_3.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 24 & plot_tibble$position >= 19] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()

pdf(file = "127_module_4.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 28 & plot_tibble$position >= 25] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()

pdf(file = "127_module_5.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 36 & plot_tibble$position >= 29] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()

pdf(file = "127_module_6.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 40 & plot_tibble$position >= 37] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()

pdf(file = "127_module_7.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 52 & plot_tibble$position >= 41] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()

pdf(file = "127_module_8.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 57 & plot_tibble$position >= 53] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()

pdf(file = "127_module_9.pdf",
    width = 4,
    height = 4)
pie(x = table(plot_tibble$psi[plot_tibble$position <= 63 & plot_tibble$position >= 58] > 49.1), init.angle = 90, col = my_colours[c(2,11)], border = NA, labels = NA)
dev.off()
#







pdf(file = "127_hist_module_1.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 9 & plot_tibble$position >= 1,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()

pdf(file = "127_hist_module_2.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 18 & plot_tibble$position >= 10,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()

pdf(file = "127_hist_module_3.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 24 & plot_tibble$position >= 19,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()

pdf(file = "127_hist_module_4.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 28 & plot_tibble$position >= 25,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()

pdf(file = "127_hist_module_5.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 36 & plot_tibble$position >= 29,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()

pdf(file = "127_hist_module_6.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 40 & plot_tibble$position >= 37,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()

pdf(file = "127_hist_module_7.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 52 & plot_tibble$position >= 41,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()

pdf(file = "127_hist_module_8.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 57 & plot_tibble$position >= 53,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()

pdf(file = "127_hist_module_9.pdf",
    width = 4,
    height = 4)
ggplot(data = plot_tibble[plot_tibble$position <= 63 & plot_tibble$position >= 58,],
       mapping = aes(x = psi)) +
  geom_vline(xintercept = 49.1, lty = 2, colour = "gray50") +
  geom_density(size = 1) +
  xlim(0,100) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)
dev.off()




















psi_vector <- c()
position_vector <- c()
k_vector <- c()

# sliding windows size 1-10
# each window deletion is defined by all the 1nt deletions that fall within the window
for (k in 1:10) {
  
  # each position
  for (p in 1:(63-k+1)) {
    
    window_start <- p
    window_end <- p+k-1
    
    temp_psi_vector <- c()
    
    for (q in window_start:window_end) {
      this_id <- paste("Del_k1_",
                       q,
                       "-",
                       q,
                       sep = "",
                       collapse = "")
      
      idx <- grep(pattern = this_id, x = indels_library$id)
      
      temp_psi_vector <- c(temp_psi_vector,
                           indels_library$psi[idx])
    }
    
    
    psi_vector <- c(psi_vector,
                    mean(temp_psi_vector))
    position_vector <- c(position_vector,
                         p + (k/2) - 0.5)
    k_vector <- c(k_vector,
                  k)
    
  }
  
}





plot_tibble <- tibble(psi = psi_vector,
                      k = k_vector,
                      position = position_vector)


ggplot(data = plot_tibble,
       mapping = aes(x = position,
                     y = psi,
                     colour = as.factor(k))) +
  geom_hline(yintercept = 49.1,
             colour = "gray80",
             size = 2) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank())































psi_vector <- c()
position_vector <- c()
k_vector <- c()

# sliding windows size 1-10
# each window deletion is defined by all the xnt deletions that fall within the window
for (k in 1:10) {
  
  # each position
  for (p in 1:(63-k+1)) {
    
    window_start <- p
    window_end <- p+k-1
    
    temp_psi_vector <- c()
    
    for (deletion_length in 1:k) {
      
      for (q in window_start:window_end) {
        this_id <- paste("Del_k",
                         deletion_length,
                         "_",
                         q,
                         "-",
                         q+deletion_length-1,
                         sep = "",
                         collapse = "")
        
        idx <- grep(pattern = this_id, x = indels_library$id)
        
        temp_psi_vector <- c(temp_psi_vector,
                             indels_library$psi[idx])
      }
      
      
    }
    
    
    
    psi_vector <- c(psi_vector,
                    mean(temp_psi_vector))
    position_vector <- c(position_vector,
                         p + (k/2) - 0.5)
    k_vector <- c(k_vector,
                  k)
    
  }
  
}





plot_tibble <- tibble(psi = psi_vector,
                      k = k_vector,
                      position = position_vector)


ggplot(data = plot_tibble,
       mapping = aes(x = position,
                     y = psi,
                     colour = as.factor(k))) +
  geom_hline(yintercept = 49.1,
             colour = "gray80",
             size = 2) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank())

