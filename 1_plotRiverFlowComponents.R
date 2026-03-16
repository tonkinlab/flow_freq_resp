##  Library packages
librarian::shelf(tidyverse, sf, ggspatial, patchwork, cowplot, lubridate, quiet = TRUE)

theme_set(theme_bw())
theme_update(text = element_text(size = 8),
             panel.grid = element_blank())


dat <- readRDS("Data/Data_models.rds")
env <- dat$EnvData %>% 
  mutate(Date = ymd(Date))

NRWQNSites <- st_read("./Data/SitesMetadata/NRWQNSiteMetadata.shp") %>% 
  filter(LocID %in% env$Site)

NRWQNSites$Latitude <- as.numeric(NRWQNSites$Latitude)

## NZ boundaries
NZReg <- st_read("./Data/reg_council_boundaries/meshblock-2020-generalised.shp")
NZReg <- NZReg[NZReg$LANDWATER_ == "Mainland",]
NZBound <- st_union(NZReg)
NZBound <- st_transform(NZBound, st_crs(NRWQNSites))

## NRWQN site data
metadata <- read.csv("./Data/SitesMetadata/NRWQNSitesMetadata.csv") %>%
  rename(site = LocID)

## Flow data
flow <- read.csv("./Data/Flow/FlowData_2022_decomp.csv") %>% 
  mutate(
    Date = ymd(Date),
    Month = month(Date,label = TRUE),
    Year = year(Date),
    Slow = slow - mean)


FlowPlot <- flow %>% 
  filter(
    site %in% env$Site,
    site != "TK3") %>%
  left_join(metadata, by = "site")


FlowPanelData <- FlowPlot %>%
  dplyr::select(Date, site, Seasonal = seas, Fast = fast, Slow, Northing) %>% 
  gather(
    key = "Component", 
    value = "Value", 
    -Date, -site, -Northing)

(flowpanels <- FlowPanelData %>%
  ggplot(aes(x = Date,
             y = Value)) +
  geom_line(aes(group = site,
                colour = Northing),
                 linewidth = 0.1) +
  facet_wrap(
    ~ Component,
    nrow = 3,
    scales = 'free_y') +
  theme(legend.position = "none", 
        strip.background = element_blank(),
        panel.border = element_blank(),
        strip.text = element_text(face = 2),
        text = element_text(size = 8),
        panel.grid.major.x = element_line(colour = 'grey80'),
        panel.grid.minor.x = element_line(colour = 'grey90'),
        axis.ticks.y = element_line()) +
  xlab("Date") +
  ylab("Log-flow") + 
  scale_colour_viridis_c() +
  scale_x_date(date_breaks = '5 year'))

##  Plotting the map
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


# VARIANCE PART
var_data <- flow %>% 
  filter(site %in% env$Site) %>%
  group_by(site) %>% 
  summarise(mean = mean(mean),
            var_seas = var(seas)/(var(seas)+var(fast)+ var(slow)),
            var_fast = var(fast)/(var(seas)+var(fast)+ var(slow)), 
            var_slow = var(slow)/(var(seas)+var(fast)+ var(slow))) %>%
  na.omit() 

mean_var_seas <- mean(var_data$var_seas)
mean_var_fast <- mean(var_data$var_fast)
mean_var_slow <- mean(var_data$var_slow)

#### FIGURE FLOW VAR 
(var_part_plot <- flow %>% 
  inner_join(metadata) %>%
  filter(site %in% env$Site) %>%
  group_by(site, Latitude) %>% 
  summarise(mean = mean(mean),
            Seas = var(seas) / (var(seas) + var(fast) + var(slow)),
            Fast = var(fast) / (var(seas) + var(fast) + var(slow)), 
            Slow = var(slow) / (var(seas) + var(fast) + var(slow))) %>%
  rename( !!paste0('Seasonal (mean = ', round(mean_var_seas,2)*100, '% )') := Seas,
          !!paste0('Fast (mean = ', round(mean_var_fast,2)*100, '% )') := Fast,
          !!paste0('Slow (mean = ', round(mean_var_slow,2)*100, '% )') := Slow) %>%
  pivot_longer(cols = -c(1:3)) %>%
  na.omit() %>%
  ggplot(aes(x = fct_reorder(site, Latitude), y = value, group = name)) +
  theme(panel.grid = element_blank(),
        text = element_text(size = 8),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.5),
        legend.position = 'inside',
        legend.position.inside = c(0.105,0.84),
        legend.key.size = unit(0.3, units = 'cm'),
        legend.background = element_rect(colour = 'grey20', linewidth = 0.2),
        legend.direction = 'vertical',
        legend.title = element_blank()) +
  geom_col(col = 'grey20', linewidth = 0.2, width = 1, aes(alpha = name, fill = Latitude)) +
  scale_y_continuous('log-flow variance (%)', expand = c(0,0), labels = scales::percent_format()) +
  guides(fill = 'none') +
  xlab('site') +
  scale_fill_viridis_c()) 


plot_grid(
  plot_grid(
    NRWQNMap,
    flowpanels,
    ncol = 2,
    rel_widths = c(2,3),
    labels = c('A', 'B')),
  var_part_plot, ncol = 1, 
  labels = c('','C'), 
  rel_heights = c(2,1.2)
  )

ggsave("Fig2FlowComponents2.pdf",
       path = "./Figs",
       device = "pdf",
       height = 16,
       width = 18,
       units = "cm",
       dpi = 300)

ggsave("Fig2FlowComponents2.png",
       path = "./Figs",
       height = 16,
       width = 18,
       units = "cm",
       dpi = 300)
