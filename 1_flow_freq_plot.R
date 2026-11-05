library(tidyverse)

imputed_flow_14 <- read.csv('Data/Flow/FlowData_2022_imputed.csv') %>%
  mutate(Date = as.Date(Date))

data <- readRDS("./Data/Data_models.RDS")
imputed_flow_14 <- imputed_flow_14 %>% filter(site %in% unique(data$Ymod$Site))

# settings 
var_name <- "Flow_input" # log flow including the seasonal cycle
n_harm <- 6  # seasonal harmonics (365, 182, 122, 91 d)
min_cycles <- 2  # drop periods with fewer than this many cycles in the record
band_cuts  <- c(30, 90, 180, 365, 2 * 365, 3 * 365)
seas_at <- 365
open_len <- 365  # length (days) drawn for the open-ended last band in the scale bar


site_analysis <- function(df, var, n_harm, min_cycles) {
  df <- df[order(df$Date), ]
  full <- tibble(Date = seq(min(df$Date), max(df$Date), by = "day")) %>%
    left_join(df[, c("Date", var)], by = "Date")

  x <- as.numeric(na.contiguous(ts(full[[var]], frequency = 1)))
  n <- length(x)
  t <- seq_len(n)
  
  # seasonal harmonics by regression: intercept + sin/cos for k = 1..n_harm (period 365.25/k)
  Xf <- cbind(1, do.call(cbind, lapply(seq_len(n_harm), function(k)
    cbind(sin(2 * pi * k * t / 365.25), cos(2 * pi * k * t / 365.25)))))
  resid      <- lm.fit(Xf, x)$residuals
  seas_share <- 1 - var(resid) / var(x)
  
  periodogram <- function(y) spec.pgram(ts(y, frequency = 1), 
                                        taper = 0.1, detrend = TRUE, plot = FALSE)
  sp_sep  <- periodogram(resid)
  sp_full <- periodogram(x)          # same length, so same frequencies as sp_sep
  
  tibble(period = 1 / sp_full$freq, share_sep = sp_sep$spec, share_full = sp_full$spec,
         seas_share = seas_share, n_days = n) %>%
    filter(period <= n / min_cycles) %>%
    mutate(share_sep  = share_sep / sum(share_sep) * (1 - seas_share),
           share_full = share_full / sum(share_full))
}


spec_df <- imputed_flow_14 %>%
  group_split(site) %>%
  set_names(map_chr(., ~ as.character(.x$site[1]))) %>%
  map_dfr(site_analysis, var = var_name, n_harm = n_harm, 
          min_cycles = min_cycles, .id = "site")

# ---- variance budget --------------------------------------------------
fmt_days <- function(d) ifelse(d < 365, paste0(d, " d"), paste0(round(d / 365, 1), " yr"))

band_levels <- c(sprintf("< %s", fmt_days(band_cuts[1])),
                 sprintf("%s - %s", fmt_days(head(band_cuts, -1)), fmt_days(tail(band_cuts, -1))),
                 sprintf("> %s", fmt_days(tail(band_cuts, 1))))
seas_level  <- "Seasonal"

# stack order, top to bottom: fast bands ... slow bands (ggplot stacks the first level on top)
n_above     <- sum(c(band_cuts, Inf) <= seas_at)     # bands above the seasonal bar
comp_levels <- append(band_levels, seas_level, after = n_above)

by_band <- spec_df %>%
  mutate(component = cut(period, c(0, band_cuts, Inf), labels = band_levels)) %>%
  group_by(site, component) %>%
  summarise(sep = sum(share_sep), full = sum(share_full), .groups = "drop") %>%
  tidyr::complete(site, component, fill = list(sep = 0, full = 0))

budget_sep <- by_band %>%
  transmute(site, component, share = sep) %>%
  bind_rows(distinct(spec_df, site, share = seas_share) %>% mutate(component = seas_level)) %>%
  mutate(component = factor(component, levels = rev(comp_levels)))

budget_full <- by_band %>%
  transmute(site, component, share = full) %>%
  mutate(component = factor(component, levels = rev(band_levels)))   # same stacking as above

comp_cols <- setNames(c(viridis::viridis(length(band_levels)), "grey"),
                      c(band_levels, seas_level))

x_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_budget_sep <- ggplot(budget_sep, aes(site, share, fill = component)) +
  geom_col() +
  scale_fill_manual(values = comp_cols, name = NULL) +
  labs(x = NULL, y = "Share of total variance",
       title = "Variance budget per site (seasonal cycle separated)",
       subtitle = paste0("Grey: first ", n_harm, " seasonal harmonics.  ",
                         "Colours: non-seasonal variance by period band.")) +
  theme_bw() + x_theme + theme(legend.position = "none")   # the scale bar below is the key

p_budget_full <- ggplot(budget_full, aes(site, share, fill = component)) +
  geom_col() +
  scale_fill_manual(values = comp_cols[band_levels], name = NULL) +
  labs(x = NULL, y = "Share of total variance",
       title = "Variance budget per site (seasonal cycle left in)",
       subtitle = "Seasonal variance is included in the period bands.") +
  theme_bw() + x_theme + theme(legend.position = "none")

#  temporal scale bar (replaces the legend) 
scale <- tibble(component = band_levels,
                   lower = c(0, band_cuts),
                   upper = c(band_cuts, max(band_cuts) + open_len),
                   open  = c(rep(FALSE, length(band_cuts)), TRUE)) %>%
  mutate(mid = (lower + upper) / 2,
         label = ifelse(open, paste0(component, " (open-ended)"), component)) %>%
  ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = lower, ymax = upper, fill = component, alpha = open),
              colour = "white") +
    geom_text(aes(x = 1.15, y = mid, label = label), hjust = 0, size = 3) +
    scale_fill_manual(values = comp_cols, guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.45), guide = "none") +
    scale_y_continuous(breaks = c(0, band_cuts), labels = c("0", fmt_days(band_cuts))) +
    scale_x_continuous(limits = c(0, 3.4)) +
    labs(x = NULL, y = "Period (segment length to scale)") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), panel.grid = element_blank(),
          axis.ticks.y = element_line())



cowplot::plot_grid(p_budget_full, scale, ncol = 2, rel_widths = c(4,1))
cowplot::plot_grid(p_budget_sep, scale, ncol = 2, rel_widths = c(4,1))
