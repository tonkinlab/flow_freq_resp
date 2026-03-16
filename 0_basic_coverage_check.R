library(ggplot2)
library(tidyverse)
library(lubridate)
library(cowplot)

# Plotting preferences
theme_set(theme_bw())
theme_update(text = element_text(size = 8),
             panel.grid = element_blank())


data = readRDS("Data/Data_models.RDS")
Y = data$Ymod
coord <- read.csv("./Data/SitesMetadata/NRWQNSitesMetadata.csv") %>% rename(Site = LocID)
site_pairs <- read.csv2("./Data/SitesMetadata/site_pairs.csv") %>% pivot_longer(-Pair, values_to = 'Site', names_to = 'Impacted')


#### N sites, samples, species ####
length(unique(data$Spt$Site)) # N sites
nrow(data$Ymod) # N samples
length(names(data$Ymod[,-c(1,2)])) # N species

#### Site characteristics ####
range(coord$Latitude) # Range of latitudes
range(data$EnvData$Order) # Range of river orders
range(coord$SiteElev) # Range of site elevations
mean(coord$SiteElev) 
median(unique(data$EnvData[,c('Site','Order')])$Order) #  mean river order

#### ORDER AND MEAN FLOW ARE HIGHLY CORRELATED ####
inner_join(data$EnvData,coord) %>% 
  ggplot(aes(y = meanFlow, x = Order)) + 
  geom_point(size = 2, shape = 21, aes(fill = Latitude)) + 
  geom_smooth(formula = y ~ x, se = F, col = 'coral') + 
  scale_fill_viridis_c('Latitude (ºN)') + 
  ylab('mean log-flow') + 
  xlab('stream order')

ggsave('Figs/correlation_order_flow.png', dpi = 600, width = 12, height = 10, units = 'cm' )
ggsave('Figs/correlation_order_flow.pdf', dpi = 600, width = 12, height = 10, units = 'cm' )

#### N unique streams #### 
N_years <- data$EnvData %>% group_by(Site) %>% summarise(n_year = length(unique(year(Date))))
hist(N_years$n_year, xlab = 'years', ylab = 'sites')


#### Seasons #### 
match_seasons <- c(rep('summer',3), rep('autumn', 4), rep('winter', 4), rep('spring', 4), 'summer')
N_season <- data$EnvData %>% mutate(season = match_seasons[month(Date)]) 
season_prop <- xtabs(~N_season$season)/nrow(N_season) * 100
paste0(names(season_prop),': ', round(season_prop, 3), '%')

#### TEMPORAL COVERAGE ####
plot_temp_cover <- inner_join(data$EnvData,coord) %>% 
  inner_join(N_years) %>%
  mutate(Site = fct_reorder(Site, Latitude),
         season =  match_seasons[month(Date)]) %>%
  ggplot(aes(y = Site, x = year(Date))) + 
  geom_tile(aes(fill = Latitude), colour = 'grey90', linewidth = 0.5) + 
  scale_x_continuous('Year',expand = c(0,0), breaks = seq(1980, 2020, 5)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(),
        panel.background = element_rect(fill = 'grey90'),
        panel.border = element_blank()) + 
  scale_fill_viridis_c('Latitude (ºN)') 

plot_n_years <- inner_join(N_years, coord) %>%
  mutate(Site = fct_reorder(Site, Latitude)) %>% 
  ggplot(aes(y = Site, x = n_year)) + 
  geom_col(aes(fill = Latitude), width = 0.7) + 
  scale_fill_viridis_c('Latitude (ºN)') + 
  xlab('sampled years') +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.title.y = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank())

plot_grid(plot_temp_cover, plot_n_years, rel_widths = c(4,1.5))

ggsave('Figs/year_coverage.png', dpi = 600, width = 18, height = 12, units = 'cm' )
ggsave('Figs/year_coverage.pdf', dpi = 600, width = 18, height = 12, units = 'cm' )

##### BASELINE VS. IMPACTED ####
inner_join(data$EnvData, site_pairs) %>% group_by(region = substr(Site, 1,2), Impacted) %>%
  summarise(SUD_N_min = min(SUD_N),
            SUD_N_max = max(SUD_N),
            SUD_N = mean(SUD_N)) %>%
  ggplot(aes(x = Impacted, y = SUD_N + 1)) + 
  theme(legend.position = 'none') +
  geom_line(aes(group = region), position = position_dodge(0.3), col = 'grey90') +
  geom_linerange(aes(ymin = SUD_N_min + 1, ymax = SUD_N_max + 1, x = Impacted, group = region, col = Impacted), 
                 position = position_dodge(0.3)) +
  geom_point(aes(col = Impacted, group = region), position = position_dodge(0.3), shape = 21, fill = 'white') + 
  scale_x_discrete(expand = c(0.2,0.2)) +
  scale_y_log10(breaks = c(1,2,11,101, 1001), labels = c(0,10^(0:3))) +
  scale_colour_manual(values = c('steelblue', 'coral')) +
  ylab('TNER') +
  xlab('site classification')
  
  
ggsave('Figs/site_pairs_TNER.pdf', dpi = 600, width = 12, height = 12, units = 'cm' )
ggsave('Figs/site_pairs_TNER.png', dpi = 600, width = 12, height = 12, units = 'cm' )
