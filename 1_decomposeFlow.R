library(tidyverse)
library(zoo)
library(lubridate)
library(imputeTS)
library(mar1s)
library(readr)
library(sf)

# function for convolution
half_gaussian_kernel <- function(sigma, size) {
  x <- seq(0, size, length.out = size + 1)
  kernel <- exp(-0.5 * (x / sigma)^2)
  kernel <- kernel / sum(kernel)  # Normalize the kernel
  return(kernel)
}

# flow data
dat <- read.csv('Data/Flow/FlowData_2022.csv')
dat$Date <- dmy(dat$Date)
 
# create list of dfs for each site
tslist <- split(dat, f = dat$site)

# Plot kernel weights 

# FLOW DECOMPOSITION 
tslist_imputed <- lapply(tslist, function(x){
  x <- x[order(x$Date), c('Date', 'Flow', 'site')]
  x_ts <- log(ts(x$Flow, frequency = 365, start = min(x$Date)))
  x_filled <- imputeTS::na_seadec(x_ts, algorithm = 'mean')
  x_mean <- mean(x_filled)
  x_seas <- rep(seasonal.smooth(x_filled - x_mean, basis = create.fourier.basis(nbasis = 10)),
                  length.out = length(x_ts))
  x_input  <- imputeTS::na_interpolation(x_ts - x_seas, maxgap = 14) + x_seas
  data.frame(x, Flow_imputed = as.numeric(x_filled), seas = as.numeric(x_seas),
             Flow_input = as.numeric(x_input))     # no trimming
})

imp <- do.call(rbind, tslist_imputed)
write.csv(imp, 'Data/Flow/FlowData_2022_imputed.csv', row.names = FALSE)
