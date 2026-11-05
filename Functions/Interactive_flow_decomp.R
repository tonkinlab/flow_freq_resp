
half_gaussian_kernel <- function(sigma, size) {
  x <- seq(0, size, length.out = size + 1)
  kernel <- exp(-0.5 * (x / sigma)^2)
  kernel / sum(kernel)
}

decompose_flow <- function(imp, kernels, seasonal_comp = TRUE, site, dates, scale_flow = FALSE,
                           slow_window = NULL, fast_window = NULL) {
  stopifnot(is.logical(seasonal_comp), length(seasonal_comp) == 1, length(site) == 1,
            length(kernels) == 2, is.logical(scale_flow), length(scale_flow) == 1)
  if (!is.null(names(kernels))) {
    if (!setequal(names(kernels), c("slow", "fast")))
      stop("named kernels must be called 'slow' and 'fast', not: ", paste(names(kernels), collapse = ", "))
    slow_size <- kernels[["slow"]]
    fast_size <- kernels[["fast"]]
  } else {                                        # unnamed: c(slow, fast)
    slow_size <- kernels[1]
    fast_size <- kernels[2]
  }
  if (slow_size <= fast_size)
    warning("slow kernel (", slow_size, ") is not larger than fast kernel (", fast_size,
            "); unnamed kernels are c(slow, fast)", call. = FALSE)

  if (is.null(slow_window)) slow_window <- slow_size
  if (is.null(fast_window)) fast_window <- fast_size
  stopifnot(slow_window >= slow_size, fast_window >= fast_size)
  dates <- as.Date(dates)
  
  needed <- c("site", "Date", "Flow_imputed", "seas")
  if (!all(needed %in% names(imp)))
    stop("imp is missing the column(s): ", paste(setdiff(needed, names(imp)), collapse = ", "),
         ". It needs: ", paste(needed, collapse = ", "), " (see prepare_imputed_flow.R)")
  
  # this site's series (Date is converted here, so imp may come from read.csv with Date as text)
  d <- imp[as.character(imp$site) == as.character(site), ]
  d$Date <- as.Date(d$Date)
  d <- d[order(d$Date), ]
  n <- nrow(d)
  if (n == 0) stop("site '", site, "' not found in imp")
  if (!is.numeric(d$Flow_imputed) || !is.numeric(d$seas))
    stop("Flow_imputed and seas must be numeric, not ", class(d$Flow_imputed)[1], " / ", class(d$seas)[1])
  if (any(diff(as.numeric(d$Date)) != 1)) stop("site '", site, "': days are not consecutive")
  burnin <- slow_window + fast_window
  if (n <= burnin)
    stop("site '", site, "': record (", n, " days) is not longer than slow + fast (", burnin, ")")
  
  x_mean <- mean(d$Flow_imputed, na.rm = TRUE)
  flow_sd <- if (scale_flow) sd(d$Flow_imputed, na.rm = TRUE) else 1
  if (!is.finite(flow_sd) || flow_sd <= 0)
    stop("site '", site, "': cannot scale, the SD of log flow is ", flow_sd)
  seas   <- if (seasonal_comp) d$seas else rep(0, n)
  x0     <- d$Flow_imputed - seas                 # series that gets convolved
  
  smooth <- function(x, size, window)               # sigma = size / 2, weights for lags 0..window
    as.numeric(stats::filter(x, half_gaussian_kernel(size / 2, window), method = "convolution", sides = 1))
  slow <- smooth(x0, slow_size, slow_window)                   # NA for the first slow_size days
  fast <- smooth(x0 - slow, fast_size, fast_window)            # NA for the first slow_size + fast_size days
  ok  <- !is.na(slow) & !is.na(fast)             # slow and fast are always NA together
  slow[!ok] <- NA_real_
  fast[!ok] <- NA_real_
  
  if (scale_flow) {
    slow <- x_mean + (slow - x_mean) / flow_sd    # keeps "slow - mean" = scaled deviation
    fast <- fast / flow_sd
    seas <- seas / flow_sd
  }
  
  m   <- match(dates, d$Date)                     # NA for dates outside the record
  out <- data.frame(site = as.character(site), Date = dates, mean = x_mean,
                    seas = if (seasonal_comp) seas[m] else NA_real_,
                    slow = slow[m], fast = fast[m], flow_sd = flow_sd)
  if ("Flow_input" %in% names(d)) out$Flow_input <- d$Flow_input[m]
  
  n_na <- sum(is.na(out$slow))
  if (n_na > 0)
    warning("site '", site, "', kernels ", slow_size, "/", fast_size, ": ", n_na, " of ",
            nrow(out), " dates have NA components (first ", burnin, " days of the record, ",
            "outside the record, or an NA in the imputed series within the last ", burnin,
            " days)", call. = FALSE)
  out
}



