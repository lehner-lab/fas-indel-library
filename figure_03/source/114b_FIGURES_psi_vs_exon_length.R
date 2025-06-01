
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



indels_library_deletions <- indels_library %>%
  filter(nchar(sequence)<=63)



exon_length <- nchar(indels_library_deletions$sequence)
exon_psi <- indels_library_deletions$psi


# This code relies on the rollapply function from the "zoo" package.  My thanks goes to Achim Zeileis and Gabor Grothendieck for their work on the package.
Quantile.cobs <- function(Y, X = NULL,
                          number.of.splits = NULL,
                          window.size = 20,
                          percent.of.overlap.between.two.windows = NULL,
                          the.distance.between.each.window = NULL,
                          the.quant = .95,
                          window.alignment = c("center"),
                          window.function = function(x) {quantile(x, the.quant)},
                          # If you wish to use this with a running average instead of a running quantile, you could simply use:
                          # window.function = mean,
                          ...)
{
  # input: Y and X, and smothing parameters
  # output: new y and x
  # Extra parameter "..." goes to the loess
  # window.size ==  the number of observation in the window (not the window length!)
  # "number.of.splits" will override "window.size"
  # let's compute the window.size:
  if(!is.null(number.of.splits)) {window.size <- ceiling(length(Y)/number.of.splits)}
  # If the.distance.between.each.window is not specified, let's make the distances fully distinct
  if(is.null(the.distance.between.each.window)) {the.distance.between.each.window <- window.size}
  # If percent.of.overlap.between.windows is not null, it will override the.distance.between.each.window
  if(!is.null(percent.of.overlap.between.two.windows))
  {
    the.distance.between.each.window <- window.size * (1-percent.of.overlap.between.two.windows)
  }
  # loading zoo
  if(!require(zoo))
  {
    print("zoo is not installed - please install it.")
    install.packages("zoo")
  }
  if(is.null(X)) {X <- index(Y)} # if we don't have any X, then Y must be ordered, in which case, we can use the indexes of Y as X.
  # creating our new X and Y
  zoo.Y <- zoo(x = Y, order.by = X)
  #zoo.X <- attributes(zoo.Y)$index
  new.Y <- rollapply(zoo.Y, width = window.size,
                     FUN = window.function,
                     by = the.distance.between.each.window,
                     align = window.alignment)
  new.X <- attributes(new.Y)$index
  new.Y.cobs <- cobs::cobs(x = new.X,
                           y = as.numeric(new.Y),
                           method = "quantile",
                           constraint="increase",
                           nknots = 10, degree = 2)$fitted
  return(list(y = new.Y, x = new.X, y.cobs = new.Y.cobs))
}






colour_function <- colorRampPalette(rev(colorspace::divergingx_hcl(n = 11, palette = "Spectral")[2:10]))
my_colours <- colour_function(11)





pdf(file = "114b_psi_vs_length_plot.pdf",
    height = 8.5,
    width = 8,
    useDingbats = F)


plot(NULL,
     las = 1,
     ylab = "percent spliced in",
     xlab = "exon length",
     main = "estimating the minimal exon length",
     xlim = c(3,63),
     ylim = c(0,100))


point_colours <- sapply(X = indels_library_deletions$id,
                        FUN = function(x){
                          this_colour <- rgb(col2rgb("gray30")[1]/255,
                                             col2rgb("gray30")[2]/255,
                                             col2rgb("gray30")[3]/255,
                                             1)
                          
                          microexon <- any(sapply(X = c("Del_k44_15-58",
                                                        "Del_k45_17-61",
                                                        "Del_k44_1-44",
                                                        "Del_k42_1-42",
                                                        "Del_k41_1-41",
                                                        "Del_k40_1-40"),
                                                  FUN = function(y){
                                                    grepl(pattern = y,
                                                          x = x)
                                                  }))
                          
                          if (microexon) {
                            this_colour <- my_colours[10]
                          }
                          
                          this_colour
                          
                        })

point_cex <- sapply(X = indels_library_deletions$id,
                    FUN = function(x){
                      this_cex <- 0.5
                      
                      microexon <- any(sapply(X = c("Del_k44_15-58",
                                                    "Del_k45_17-61",
                                                    "Del_k44_1-44",
                                                    "Del_k42_1-42",
                                                    "Del_k41_1-41",
                                                    "Del_k40_1-40"),
                                              FUN = function(y){
                                                grepl(pattern = y,
                                                      x = x)
                                              }))
                      
                      if (microexon) {
                        this_cex <- 1
                      }
                      
                      this_cex
                      
                    })

points(x = exon_length,
       y = exon_psi,
       pch = 19,
       cex = point_cex,
       col = point_colours)

abline(v = 27)
abline(h=3)

View(indels_library_deletions %>%
       filter(nchar(sequence) <=27) %>%
       filter(psi > 3))
for (i in 1:5) {
  
  q <- seq(0.25,0.65,0.1)[i]
  
  cobs_model_1 <- Quantile.cobs(Y = exon_psi,
                                X = exon_length,
                                window.size = 20,
                                the.quant = q)
  cobs_model_2 <- Quantile.cobs(Y = exon_psi,
                                X = exon_length,
                                window.size = 20,
                                the.quant = q+0.1)
  
  polygon(x = c(cobs_model_1$x,
                rev(cobs_model_2$x)),
          y = c(cobs_model_1$y.cobs,
                rev(cobs_model_2$y.cobs)),
          border = rgb(col2rgb(my_colours[2+i])[1]/255,
                       col2rgb(my_colours[2+i])[2]/255,
                       col2rgb(my_colours[2+i])[3]/255,
                       0.75),
          col = rgb(col2rgb(my_colours[2+i])[1]/255,
                    col2rgb(my_colours[2+i])[2]/255,
                    col2rgb(my_colours[2+i])[3]/255,
                    0.75)
  )
}

cobs_model <- Quantile.cobs(Y = exon_psi,
                            X = exon_length,
                            window.size = 20,
                            the.quant = 0.5)

lines(x = cobs_model$x,
      y = cobs_model$y.cobs,
      lwd = 3,
      col = "black")
dev.off()



























pdf(file = "114b_psi_vs_length_plot_IQR.pdf",
    height = 8.5,
    width = 8,
    useDingbats = F)


plot(NULL,
     las = 1,
     ylab = "percent spliced in",
     xlab = "exon length",
     main = "estimating the minimal exon length",
     xlim = c(3,63),
     ylim = c(0,100))


point_colours <- sapply(X = indels_library_deletions$id,
                        FUN = function(x){
                          this_colour <- rgb(col2rgb("gray30")[1]/255,
                                             col2rgb("gray30")[2]/255,
                                             col2rgb("gray30")[3]/255,
                                             0.75)
                          
                          microexon <- any(sapply(X = c("Del_k44_15-58",
                                                        "Del_k45_17-61",
                                                        "Del_k44_1-44",
                                                        "Del_k42_1-42",
                                                        "Del_k41_1-41",
                                                        "Del_k40_1-40"),
                                                  FUN = function(y){
                                                    grepl(pattern = y,
                                                          x = x)
                                                  }))
                          
                          if (microexon) {
                            this_colour <- my_colours[10]
                          }
                          
                          this_colour
                          
                        })

point_cex <- sapply(X = indels_library_deletions$id,
                    FUN = function(x){
                      this_cex <- 0.5
                      
                      microexon <- any(sapply(X = c("Del_k44_15-58",
                                                    "Del_k45_17-61",
                                                    "Del_k44_1-44",
                                                    "Del_k42_1-42",
                                                    "Del_k41_1-41",
                                                    "Del_k40_1-40"),
                                              FUN = function(y){
                                                grepl(pattern = y,
                                                      x = x)
                                              }))
                      
                      if (microexon) {
                        this_cex <- 1
                      }
                      
                      this_cex
                      
                    })

points(x = exon_length,
       y = exon_psi,
       pch = 19,
       cex = point_cex,
       col = point_colours)


cobs_model_1 <- Quantile.cobs(Y = exon_psi,
                              X = exon_length,
                              window.size = 20,
                              the.quant = 0.25)
cobs_model_2 <- Quantile.cobs(Y = exon_psi,
                              X = exon_length,
                              window.size = 20,
                              the.quant = 0.75)

polygon(x = c(cobs_model_1$x,
              rev(cobs_model_2$x)),
        y = c(cobs_model_1$y.cobs,
              rev(cobs_model_2$y.cobs)),
        border = rgb(col2rgb(my_colours[3])[1]/255,
                     col2rgb(my_colours[3])[2]/255,
                     col2rgb(my_colours[3])[3]/255,
                     0.75),
        col = rgb(col2rgb(my_colours[3])[1]/255,
                  col2rgb(my_colours[3])[2]/255,
                  col2rgb(my_colours[3])[3]/255,
                  0.75))
cobs_model <- Quantile.cobs(Y = exon_psi,
                            X = exon_length,
                            window.size = 20,
                            the.quant = 0.5)

lines(x = cobs_model$x,
      y = cobs_model$y.cobs,
      lwd = 3,
      col = "black")
dev.off()

write_delim(x = tibble(exon_length = exon_length,
                       psi = exon_psi),
            delim = "\t",
            file = "114b_figure_03a.tsv")
