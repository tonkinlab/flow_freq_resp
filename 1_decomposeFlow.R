library(tidyverse)
library(zoo)
library(lubridate)
library(imputeTS)
library(mar1s)
library(readr)
library(sf)

# convolution kernel sizes (in days, approx. 95% of density)
slow_kernel_size = 365*2 
fast_kernel_size = 90 

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
times = 1:730
slow_kernel <- half_gaussian_kernel(slow_kernel_size/2, 729)
fast_kernel <- half_gaussian_kernel(fast_kernel_size/2, 729)
plot(times, fast_kernel/max(fast_kernel), type = 'l')
lines(times, slow_kernel/max(slow_kernel),  col = 'red')
abline(v = slow_kernel_size, col = 'red')
abline(v = fast_kernel_size)

# FLOW DECOMPOSITION 
tslist_input_14 <- lapply(tslist, function(x){
  x <- x[,c('Date', 'Flow', 'site')]
  x_ts <- log(ts(x$Flow, frequency = 365, start = min(x$Date)))
  x_seas_imput <-  imputeTS::na_seadec(x_ts, algorithm = 'mean')
  x_mean <- mean(x_seas_imput)
  x_seas <- rep(seasonal.smooth(x_seas_imput - x_mean, basis = create.fourier.basis(nbasis = 10)), length.out = length(x_ts))
  # Interpolate gaps 
  x_imput <- imputeTS::na_interpolation(x_ts - x_seas, maxgap = 14) + x_seas
  # consecutive convolutions (slow,  fast)
  x_slow_2m <- c(rep(x_mean, slow_kernel_size), 
                 convolve(x_seas_imput - x_seas, rev(half_gaussian_kernel(slow_kernel_size/2, slow_kernel_size)), type = 'filter'))
  x_fast <-  c(rep(-9999, fast_kernel_size), 
               convolve(unclass(x_seas_imput - x_seas - x_slow_2m), rev(half_gaussian_kernel(fast_kernel_size/2, fast_kernel_size)), type = 'filter'))
  data.frame(x[-c(1:(slow_kernel_size + fast_kernel_size)),],
             mean = x_mean, 
             seas = unclass(x_seas)[-c(1:(slow_kernel_size + fast_kernel_size))], 
             fast = unclass(x_fast)[-c(1:(slow_kernel_size + fast_kernel_size))],
             slow = unclass(x_slow_2m)[-c(1:(slow_kernel_size + fast_kernel_size))],
             Flow_input = unclass(x_imput)[-c(1:(slow_kernel_size + fast_kernel_size))])
})

# converting to df 
imputed_flow_14 <- data.frame(do.call(rbind,tslist_input_14))
write.csv(imputed_flow_14, 'Data/Flow/FlowData_2022_decomp.csv', row.names = FALSE)
