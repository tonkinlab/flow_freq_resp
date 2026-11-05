# Var plot:
library(tidyverse)
library(sf)
library(ggspatial)
library(ggrepel)
library(cowplot)
library(lubridate)

source('Functions/Interactive_flow_decomp.R')   # half_gaussian_kernel() and decompose_flow()

theme_set(theme_bw(base_size = 8))
theme_update(panel.grid = element_blank())

imputed_flow_14 <- read.csv('Data/Flow/FlowData_2022_imputed.csv') %>%
  mutate(Date = as.Date(Date))

data <- readRDS("./Data/Data_models.RDS")
imputed_flow_14 <- imputed_flow_14 %>% filter(site %in% unique(data$Ymod$Site))

env <- data$EnvData %>%
  mutate(Date = ymd(Date))

## NRWQN site data
metadata <- read.csv("./Data/SitesMetadata/NRWQNSitesMetadata.csv") %>%
  rename(site = LocID)

# settings
var_name <- "Flow_input" # log flow including the seasonal cycle
n_harm <- 6  # seasonal harmonics (365, 182, 122, 91, 73, 61 d)
min_cycles <- 2  # drop periods with fewer than this many cycles in the record
band_cuts  <- c(30, 90, 180, 365, 2 * 365, 3 * 365)
seas_at <- 365
open_len <- 365  # length (days) drawn for the open-ended last band in the scale bar

# flow decomposition shown in panel B (kernel sizes in days)
kernels <- c(slow = 2 * 365, fast = 90)
seasonal_comp <- FALSE   # TRUE: seasonal cycle removed before slow / fast are computed


# ---- flow decomposition (panel B) -----------------------------------------
# decompose_flow() works on one site; here every date of every site is requested.
# The first slow + fast days of each record are the burn-in (NA) and are dropped.
flow <- map_dfr(unique(imputed_flow_14$site), function(s)
  suppressWarnings(
    decompose_flow(imputed_flow_14, kernels, seasonal_comp, s,
                   imputed_flow_14$Date[imputed_flow_14$site == s]))) %>%
  filter(!is.na(slow), !is.na(fast)) %>%
  mutate(
    Month = month(Date, label = TRUE),
    Year = year(Date),
    Slow = slow - mean)


# ---- variance partitioning (panel C) ---------------------------------------
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

# mean share of each component across sites (for the text)
budget_mean <- bind_rows(sep = budget_sep, full = budget_full, .id = "version") %>%
  group_by(version, component) %>%
  summarise(mean_share = mean(share), .groups = "drop")
print(budget_mean, n = Inf)

comp_cols <- setNames(c(viridis::viridis(length(band_levels)), "grey"),
                      c(band_levels, seas_level))

# sites ordered by latitude, as in the previous version of panel C
site_order <- metadata %>%
  filter(site %in% unique(spec_df$site)) %>%
  arrange(Latitude) %>%
  pull(site)
site_order <- c(site_order, setdiff(unique(spec_df$site), site_order))


# FOR SI
p_budget_sep <- ggplot(budget_sep, aes(factor(site, levels = site_order), share, fill = component)) +
  geom_col() +
  scale_fill_manual(values = comp_cols, name = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none',
        axis.ticks.y = element_line()) +
  xlab('site') +
  scale_y_continuous('log-flow variance (%)', expand = c(0,0), limits = c(0,1.01),
                     labels = scales::percent_format())

p_budget_full <- ggplot(budget_full, aes(factor(site, levels = site_order), share, fill = component)) +
  geom_col() +
  scale_fill_manual(values = comp_cols[band_levels], name = NULL) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none',
        axis.ticks.y = element_line()) +
  xlab('site') +
  scale_y_continuous('log-flow variance (%)', expand = c(0,0), limits = c(0,1.01),
                     labels = scales::percent_format())

#  temporal scale bar (replaces the legend)
p_scale <- tibble(component = band_levels,
                  lower = c(0, band_cuts),
                  upper = c(band_cuts, max(band_cuts) + open_len)) %>%
  mutate(mid = (lower + upper) / 2,
         label = component) %>%
  ggplot() +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = lower, ymax = upper, fill = component),
            colour = "white") +
  geom_text(aes(x = 1.15, y = mid, label = label), hjust = 0, size = 2) +
  scale_fill_manual(values = comp_cols, guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.45), guide = "none") +
  scale_y_continuous(breaks = c(0, band_cuts), labels = c("", fmt_days(band_cuts))) +   # no "0" label: it overlaps "30 d"
  scale_x_continuous(limits = c(0, 3.4)) +
  labs(x = NULL, y = "Period (segment length to scale)") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        axis.ticks.y = element_line())

plot_var_part <- cowplot::plot_grid(p_budget_full, p_scale, ncol = 2, rel_widths = c(4,1))
plot_var_part_seas <- cowplot::plot_grid(p_budget_sep, p_scale, ncol = 2, rel_widths = c(4,1))

# For SI
ggsave('Figs/plot_vp_flow_seas.pdf',plot_var_part_seas, width = 18, height = 10, units = 'cm', dpi = 600,
       bg = 'white')


# ---- map data ---------------------------------------------------------------
NRWQNSites <- st_read("./Data/SitesMetadata/NRWQNSiteMetadata.shp") %>%
  filter(LocID %in% env$Site)

NRWQNSites$Latitude <- as.numeric(NRWQNSites$Latitude)

## NZ boundaries
NZReg <- st_read("./Data/reg_council_boundaries/meshblock-2020-generalised.shp")
NZReg <- NZReg[NZReg$LANDWATER_ == "Mainland",]
NZBound <- st_union(NZReg)
NZBound <- st_transform(NZBound, st_crs(NRWQNSites))


##  Plotting the map (with site labels)
(NRWQNMap <- ggplot() +
    theme_map() +
    geom_sf(
      NZBound,
      mapping = aes(),
      fill = "grey90") +
    geom_sf(
      NRWQNSites,
      mapping  = aes(fill = abs(Northing)),
      color = "black",
      shape = 21) +
    ggrepel::geom_text_repel(
      data = NRWQNSites,
      aes(label = LocID, geometry = geometry),
      stat = "sf_coordinates",
      size = 1.8,
      segment.size = 0.15,
      min.segment.length = 0,
      box.padding = 0.15,
      max.overlaps = Inf,
      seed = 1) +
    coord_sf(crs = st_crs(NRWQNSites),
             xlim = c(1009050, 2050000),
             ylim = c(4800000, 6170000),
             datum = st_crs(NRWQNSites)) +
    ggspatial::annotation_north_arrow(location = "br",
                                      which_north = "magnetic",
                                      pad_y = unit(0.3, "in"),
                                      pad_x = unit(-0.1, "in"),
                                      style = ggspatial::north_arrow_minimal()) +
    ggspatial::annotation_scale(location = "br",
                                line_width = 0.5,
                                width_hint = 0.4,
                                bar_cols = c("black", "grey90")) +
    viridis::scale_fill_viridis(direction = 1) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(color = "Northing",
         size = "Mean Flow") +
    theme(text = element_text(size = 8),
          legend.position = ""))


#figure
fig2 <- plot_grid(
  NRWQNMap,
  plot_var_part,
  ncol = 2,
  rel_widths = c(1,2),
  labels = c('A','B')
)

ggsave("Fig2FlowComponents2.pdf",
       plot = fig2,
       path = "./Figs",
       device = "pdf",
       height = 10,
       width = 18,
       units = "cm",
       dpi = 300,
       bg = "white")   # the map and panel C have transparent backgrounds

ggsave("Fig2FlowComponents2.png",
       plot = fig2,
       path = "./Figs",
       height = 10,
       width = 18,
       units = "cm",
       dpi = 300,
       bg = "white")
