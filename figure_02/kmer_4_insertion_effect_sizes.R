
library(tidyverse)
library(scales)
library(pheatmap)
load("data/indels_library.RData")



idx_singles <- which(nchar(indels_library$sequence)==63 & indels_library$id != "WT")
density_object <- density(x = indels_library$es[idx_singles])
library_bias <- density_object$x[which(density_object$y == max(density_object$y))]
# plot(density_object)
# abline(v=library_bias)
# plot(density(indels_library$es - library_bias))
# abline(v=0)
# abline(v=1)
indels_library$es <- indels_library$es - library_bias

indels_library$psi <- (49.1/1) * exp(indels_library$es)

indels_library$psi[which(indels_library$psi > 100)] <- 100


indels_library <- separate_rows(indels_library, id, sep = ";")


k <- 4
idx_k <- which(nchar(indels_library$sequence)==63+k)



mean(abs(indels_library$psi[idx_k] - 49.1) > 10, na.rm = TRUE)
mean(abs(indels_library$psi[idx_k] - 49.1), na.rm = TRUE)
mean((indels_library$psi[idx_k] - 49.1) > 10, na.rm = TRUE)
mean((indels_library$psi[idx_k] - 49.1) < -10, na.rm = TRUE)


