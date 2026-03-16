library(ggplot2)
library(cowplot)
library(tidyverse)

theme_set(theme_bw())
theme_update(text = element_text(size = 8),
             panel.grid = element_blank())

#### ORDINATION OF ABUNDANCE ONLY MODEL  ####

# Read model file and check convergence
mAb <- readRDS('./Models/mAb.rds')
mAb$convergence

# Load species taxa groups and sensitivity scores
source('Data/taxa_tolerance.R')
source('Data/taxa_groups.R')


coord <- data.frame(gllvm::getLV(mAb))
coord$site <- rownames(coord)
coord <- read.csv("./Data/SitesMetadata/NRWQNSitesMetadata.csv") %>%
  rename(site = LocID) %>% dplyr::select(site, Latitude, Longitude) %>% 
  inner_join(data.frame(coord))

load <- data.frame(gllvm::getLoadings(mAb))
load$sp <- rownames(load)
load$abs <- sqrt(rowSums(load[,c(1,2)]^2))
load <- load[rev(order(load$abs)),][,]
load$group <- taxa_groups[load$sp]
load$tolerance <- tolerance[load$sp]

plot_site_coord <- ggplot(data.frame(coord), aes(x = LV1, y = LV2)) + 
  geom_hline(yintercept = 0, lty = 5, col = 'grey') + 
  geom_vline(xintercept = 0, lty = 5, col = 'grey') + 
  geom_point(shape = 21, aes(fill = Latitude)) +
  geom_text(aes(label = site, y = LV2 + 0.07), size = 2) +
  scale_fill_viridis_c() +
  ylim(range(coord$LV2)*1.05) + 
  xlim(range(coord$LV1)) + 
  coord_equal()

plot_sp_load <- ggplot(load, aes(x = LV1*0.6, y = LV2*0.6)) + 
  geom_point(data = coord, aes(x = LV1, y = LV2), col = 'grey') +
  geom_hline(yintercept = 0, lty = 5, col = 'grey') + 
  geom_vline(xintercept = 0, lty = 5, col = 'grey') + 
  geom_text(aes(label = sp, colour = tolerance), size = 1) +
  ylim(range(coord$LV2)*1.05) + 
  xlim(range(coord$LV1)) + 
  ylab('LV2') +
  xlab('LV1') +
  coord_equal() + 
  scale_colour_gradientn('sensitivity', colours = c('red', 'grey', 'steelblue')) 

(ord_plot <- plot_grid(plot_site_coord, plot_sp_load, align = 'hv', labels = c('A', 'B')))
ggsave('Figs/ord_plot.png', width = 18, height = 10, dpi = 600, units = 'cm')
ggsave('Figs/ord_plot.pdf', width = 18, height = 10, dpi = 600, units = 'cm')


#### ORDINATION OF MODEL WITH TRAITS ####

# Load models and check convergence 
mTr <- readRDS('./Models/mTr.rds')
mTr$convergence

# Ordination
#### PLOT ORDINATION ####
coord <- data.frame(gllvm::getLV(mTr))
coord$site <- rownames(coord)
coord <- read.csv("./Data/SitesMetadata/NRWQNSitesMetadata.csv") %>%
  rename(site = LocID) %>% dplyr::select(site, Latitude, Longitude) %>% inner_join(data.frame(coord))

load <- data.frame(gllvm::getLoadings(mTr))
load$sp <- rownames(load)
load$abs <- sqrt(rowSums(load[,c(1,2)]^2))
load <- load[rev(order(load$abs)),][,]
load$group <- taxa_groups[load$sp]
load$tolerance <- tolerance[load$sp]

plot_site_coord <- ggplot(coord, aes(x = LV1, y = LV2)) + 
  geom_hline(yintercept = 0, lty = 5, col = 'grey') + 
  geom_vline(xintercept = 0, lty = 5, col = 'grey') + 
  geom_point(shape = 21, aes(fill = Latitude)) +
  geom_text(aes(label = site, y = LV2 + 0.07), size = 2) +
  scale_fill_viridis_c() +
  coord_equal()

plot_sp_load <- ggplot(load, aes(x = LV1*0.6, y = LV2*0.6)) + 
  geom_point(data = coord, aes(x = LV1, y = LV2), col = 'grey') +
  geom_hline(yintercept = 0, lty = 5, col = 'grey') + 
  geom_vline(xintercept = 0, lty = 5, col = 'grey') + 
  geom_text(aes(label = sp, colour = tolerance), size = 1) +
  scale_fill_viridis_c() + 
  ylab('LV2') +
  xlab('LV1') +
  coord_equal() + 
  scale_colour_gradientn(colours = c('red', 'grey', 'steelblue')) 

(ord_plot_traits <- plot_grid(plot_site_coord, plot_sp_load, align = 'hv', labels = c('A', 'B')))

ggsave('Figs/ord_plot_traits.png', width = 18, height = 10, dpi = 600, units = 'cm')
ggsave('Figs/ord_plot_traits.pdf', width = 18, height = 10, dpi = 600, units = 'cm')


