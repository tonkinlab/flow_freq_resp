library(tidyverse)
library(ggplot2)
library(cowplot)

theme_set(theme_bw())
theme_update(text = element_text(size = 8),
             panel.grid = element_blank())


TrData <- readRDS('./Data/TrData_model.rds')

source('Functions/predict_newtraits.R')
source('Functions/trait_slope.R')

mTr <- readRDS('./Models/mTr.rds')
mTr$convergence

tr  <- mTr$TR                                                   # species x traits
b_f <- mTr$fourth.corner[preds_flow, , drop = FALSE]            # flow terms x traits
tr  <- tr[, colnames(b_f), drop = FALSE]                        # same trait order as b_f
b_e <- mTr$params$B[preds_flow]                                 # common (main) slopes
z_i <- mTr$params$Br[preds_flow, , drop = FALSE]                # species random slopes (flow terms x species)

# Pseudo R2 for slopes
re <- t(z_i)                                                    # species x flow terms
eta_marginal    <- sweep(tr %*% t(b_f), 2, b_e, FUN = '+')
eta_conditional <- sweep(tr %*% t(b_f) + t(z_i), 2, b_e, FUN = '+')


r2_trait <- c(slow = cor(eta_marginal[,1], eta_conditional[,1])^2, fast = cor(eta_marginal[,2], eta_conditional[,2])^2)

(r2_plot <- data.frame(var =  names(r2_trait), r2 = as.numeric(r2_trait)) %>%
    ggplot(aes(x = var, y = r2)) +
    geom_col(fill = 'grey', col = 'black', linewidth = 0.3) +
    theme(panel.border = element_blank(),
          axis.line.y = element_line(linewidth = 0.3),
          axis.ticks.x = element_blank(),
          aspect.ratio = 1.5) +
    xlab('flow component') +
    scale_y_continuous('Pseudo-R2 (slopes)', limits = c(0, 0.6), expand = c(0, 0)))


# 4th corner coefficients: s
summary(mTr)
fourth_c <- data.frame(t(mTr$fourth.corner))
fourth_c$trait <- row.names(fourth_c)
fourth_c_plot <- fourth_c %>%
  pivot_longer(cols = -trait, names_to = 'var', values_to = 'val') %>%
  mutate(
    # Keeping significance logic
    sig = case_when(
      trait == 'DissPotential_LOW2' & var == 'Mean' ~ '0.069 .',
      trait == 'CycPerYear_PLURIV2' & var == 'Mean' ~ '0.010 *',
      trait == 'aq_ADUANDLAR2' & var == 'Trend2' ~ '0.021 *',
      trait == 'aq_ADUORLAR2' & var == 'Trend2' ~ '0.002 **',
      trait == 'aq_ADUORLAR2' & var == 'Seasonal' ~ '0.034 *',
      trait == 'Mob_BURROWER2' & var == 'Seasonal' ~ '0.006 **',
      trait == 'Mob_SWIMMER2' & var == 'Residual' ~ '. 0.052 .',
      trait == 'aq_ADUANDLAR2' & var == 'l_SUD_N' ~ '. 0.074 .',
      TRUE ~ ''
    ),
    trait2 = recode(trait, 
                    Mob_SWIMMER2 = 'MOBILITY: swimmer', 
                    Mob_BURROWER2 = 'MOBILITY: burrower',
                    Mob_ATTACHED2 = 'MOBILITY: attached',
                    aq_ADUANDLAR2 = 'AQUATIC STAGE: adult and larva',
                    aq_ADUORLAR2 = 'AQUATIC STAGE: adult or larva',
                    DissPotential_LOW2 = 'DISSEMINATION: low',
                    DissPotential_HIGH2 = 'DISSEMINATION: high',
                    CycPerYear_SEMI2 = 'CYCLES: semivoltine',
                    CycPerYear_PLURIV2 = 'CYCLES: plurivoltine',
                    l_SIZE = 'SIZE'),
    var2 = recode(var, 
                  Mean = 'mean', 
                  Trend2 = 'slow',
                  Residual = 'fast',
                  Seasonal = 'seasonal',
                  l_SUD_N = 'TNER')
  ) %>%
  filter(var2 %in% c('slow', 'fast', 'seasonal')) %>%
  ggplot(aes(x = var2, y = trait2)) + 
  geom_tile(aes(fill = val), col = 'white', linewidth = 1) + 
  geom_text(aes(label = sig), size = 1.5) + 
  scale_fill_gradientn('coefficient (trait : component)', 
                       colours = c('steelblue','lightblue3','white','coral', 'coral2'), 
                       limits = c(-0.2, 0.2)) +
  theme(legend.title = element_text(hjust = 0.5, vjust = 1, size = 6),
        legend.key.height = unit(0.5, 'lines'),
        legend.key.width = unit(1, 'lines'),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'bottom',
        legend.title.position = 'top') 

out_trait <- list()

for (var in c('Trend2', 'Residual')) {
  for (trait in c('aq', 'DissPotential', 'Mob', 'l_SIZE')) {
    print(trait)
    if (substr(trait,1,2) %in% c('Mo', 'aq', 'Di')) {
      type = 'categorical'
    }else{
      type = 'numeric'
    }
    pred_trait <-
      trait_effects(mTr,
                    env_name = var,
                    trait_prefix = trait,
                    trait_type = type,
                    trait_value = NULL,
                    traits = TrData,
                    nsim = 10000)
    
    out_trait[[length(out_trait) + 1]] <- pred_trait
  }
}

out_trait_data <- data.frame(do.call(rbind, out_trait))
out_trait_data$baseline <- ifelse(out_trait_data$trait_category == 'baseline', T, F)
out_trait_data$trait_category <- ifelse(out_trait_data$trait_category == 'baseline' & out_trait_data$trait == 'aq', 'Larva and Pupa', out_trait_data$trait_category)
out_trait_data$trait_category <- ifelse(out_trait_data$trait_category == 'baseline' & out_trait_data$trait == 'DissPotential', 'medium', out_trait_data$trait_category)
out_trait_data$trait_category <- ifelse(out_trait_data$trait_category == 'baseline' & out_trait_data$trait == 'Mob', 'crawler', out_trait_data$trait_category)
out_trait_data$trait_category <- as.factor(out_trait_data$trait_category)

levels(out_trait_data$trait_category) <- c('adult and larva', 'adult or larva', 'attached', 'burrower', 'crawler', 'high','larva and pupa', 'low', 'max', 'mean', 'medium', 'min', 'swimmer')
out_trait_data$trait_category <- paste0(out_trait_data$trait_category, ifelse(out_trait_data$baseline, '*', ''))
out_trait_data$trait <- as.factor(out_trait_data$trait)

levels(out_trait_data$trait) <- c('aquatic\n stage', 'cycles\n per year', 'dissemination\n potential', 'size', 'mobility')
out_trait_data$env_var <- as.factor(out_trait_data$env_var)
levels(out_trait_data$env_var) <- c('fast', 'slow')

trait_slopes_plot <- out_trait_data %>% filter(env_var %in% c('fast', 'slow')) %>% 
  ggplot(aes(x = median, y = trait_category)) + 
  theme(strip.background = element_blank(), strip.placement = 'outside',
        axis.title = element_blank(),
        axis.text.x = element_text(size = 5),
        panel.border = element_rect(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = unit(c(0.2,0.5,0.2,0.7), "cm"),
        panel.spacing = unit(0.1, "lines")) + 
  scale_x_continuous(breaks = c(-1.5, -1, -0.5, 0, 0.5, 1)) + 
  geom_vline(xintercept = 0, lty = 5, linewidth = 0.3, colour = 'grey50') +
  geom_linerange(aes(xmin = lower_95, xmax = upper_95), linewidth = 0.4, colour = 'grey20') +
  geom_point(shape = 21, fill = 'white') +
  facet_grid(trait ~ env_var, scales = 'free', switch = 'x')

plot_grid(fourth_c_plot, trait_slopes_plot, r2_plot, rel_widths = c(1,1.2,0.6), labels = c('A', 'B', 'C'), ncol = 3) 

ggsave('Figs/fourt_plot.png', width = 18, height = 9, dpi = 600, units = 'cm')
ggsave('Figs/fourt_plot.pdf', width = 18, height = 9, dpi = 600, units = 'cm')
