library(ggplot2)
library(cowplot)
library(tidyverse)
library(GGally)

theme_set(theme_bw())
theme_update(text = element_text(size = 8),
             panel.grid = element_blank())

mAb <- readRDS('./Models/mAb.rds')
mAb$convergence

source('Data/taxa_groups.R')
source('Data/taxa_tolerance.R')

#### Coefficient plot ###
Xcoef <- data.frame(mAb$params[['Xcoef']]) %>% 
  mutate(sp = row.names(.)) %>% 
  pivot_longer(cols = -sp, names_to = 'var', values_to = 'coef')

coef_data <- data.frame(mAb$sd[['Xcoef']]) %>% 
  mutate(sp = row.names(.)) %>% 
  pivot_longer(cols = -sp, names_to = 'var', values_to = 'sd') %>%
  inner_join(Xcoef) %>%
  mutate(lower = coef - 1.96*sd,
         upper = coef + 1.96*sd, 
         sig = sign(lower) == sign(upper)) %>%
  mutate(var2 = recode(var, Mean = " mean", Trend2 = "slow", Residual = 'fast', Seasonal = 'seasonal', l_SUD_N = 'TNER'))

## Mean absolute responses ## 
coef_data %>% group_by(var) %>% summarise(mean_coef = mean(abs(coef)), 
                                          mean_coef_sd = sd(abs(coef)))

# Flow 
group_colours =   c('purple4','steelblue', 'darkolivegreen4', 'khaki3', 'gold2','darkorange', 'coral2', 'brown','grey20','grey20','grey20')

coef_data %>% 
  filter(var2 %in% c('fast', 'slow', 'seasonal')) %>%
  mutate(
    group = factor(taxa_groups[sp], levels = c('Ephemeroptera', 'Plecoptera', 'Trichoptera', 'Diptera', 'Coleoptera', 'Gastropoda','Crustacea', 'Worms', 'Arachnida', 'Hemiptera', 'Megaloptera')),
    var2_label = recode(var2, 
                        fast = "atop(bold('Fast'), 'MAE: 0.58')",
                        slow = "atop(bold('Slow'), 'MAE: 1.64')",
                        seasonal = "atop(bold('Seasonal'), 'MAE: 0.83')")
  ) %>%
  ggplot(aes(x = coef, y = sp)) +
  # Keep label_parsed to interpret the atop() and bold() commands
  facet_grid(group ~ var2_label, scales = 'free', space = 'free_y', labeller = label_parsed) + 
  theme(
    legend.position = 'bottom',
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.5),
    axis.line.x = element_line(linewidth = 0.5),
    legend.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(colour = 'grey95', linewidth = 0.2),
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_blank(),
    # element_text is safe for the S7 system!
    strip.text.x = element_text(colour = 'grey20', size = 10, face = "plain", lineheight = 1.1),
    strip.text.y = element_text(angle = 0, face = "plain", size = 8, hjust = 0),
    panel.border = element_rect(colour = NA, fill = NA),
    panel.background = element_rect(fill = 'white', colour = NA)
  ) +
  geom_vline(xintercept = 0, lty = 2, col = 'grey50', linewidth = 0.2) +
  geom_linerange(aes(xmin = lower, xmax = upper, colour = group, alpha = sig), linewidth = 0.4) +
  geom_point(aes(col = group, alpha = sig), stroke = 0.3) + 
  scale_alpha_manual('', values = c(0.2, 1)) +
  scale_y_discrete('') + 
  scale_colour_manual(values = group_colours) +
  guides(alpha = 'none', colour = 'none', fill = 'none') +
  xlab('Response (log-abundance / log-flow)')

ggsave('./Figs/coefs_flow.png', width = 18, height = 20, units = 'cm', dpi = 600)
ggsave('./Figs/coefs_flow.pdf', width = 18, height = 20, units = 'cm', dpi = 600)

# TNER and mean flow
coef_data %>% filter(!(var2 %in% c('fast', 'slow', 'seasonal'))) %>%
  mutate(group = taxa_groups[sp],
         sensitivity = tolerance[sp],
         sp = fct_reorder(sp, sensitivity))  %>%
  mutate() %>%
  ggplot(aes(x = coef, y = sp, group = group)) +
  facet_wrap(~var2, ncol = 5, scales = 'free_x') + 
  theme(legend.position = 'right',
        legend.title = element_text(hjust = 0.5),
        panel.grid.major.y = element_line(colour = 'grey95', linewidth = 0.2),
        panel.spacing = unit(-0.1, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = 'bold', colour = 'grey20', size = 12),
        panel.border = element_rect(linewidth = 0.5),
        panel.background = element_rect(fill = 'white')) +
  geom_vline(xintercept = 0, lty = 2, col = 'grey50') +
  geom_linerange(aes(xmin = lower, xmax = upper, alpha = sig), linewidth = 0.4) +
  geom_point(aes(alpha = sig, fill = as.numeric(sensitivity)), stroke = 0.3, shape = 21) + 
  scale_alpha_manual('', values = c(0.2,1)) +
  scale_y_discrete('Taxa') +
  scale_x_continuous(expand = c(0, 0.1)) + 
  scale_fill_gradientn('sensitivity', colours = c('red', 'pink', 'grey', 'skyblue', 'steelblue4'), limits = c(0,10)) +
  guides(alpha = 'none') +
  xlab('slope')

ggsave('./Figs/coefs_TNER.png', width = 15, height = 17, units = 'cm', dpi = 600)
ggsave('./Figs/coefs_TNER.pdf', width = 15, height = 17, units = 'cm', dpi = 600)

### COEF CORRELATION ###
data.frame(mAb$params$Xcoef) %>% rename(mean = 'Mean', slow = 'Trend2', seasonal = 'Seasonal', fast = 'Residual') %>%
  dplyr::select(mean, slow, seasonal, fast) %>%
  ggpairs( 
    diag = list(continuous = wrap("densityDiag", fill = 'wheat', col = 'black')), 
    lower = list(continuous = "points", shape = 21)) + 
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave('Figs/coef_cor.png', dpi = 600, width = 13, height = 12, units = 'cm')
ggsave('Figs/coef_cor.pdf', dpi = 600, width = 13, height = 12, units = 'cm')
