library(ggplot2)
library(cowplot)
library(gllvm)
library(tidyverse)

theme_set(theme_bw())
theme_update(text = element_text(size = 8),
             panel.grid = element_blank())

# Read model file and check convergence
mAb <- readRDS('./Models/mAb.rds')
mAb$convergence

# Load model matrix
mod.mat <- readRDS('./Data/model_matrix_env.rds')
StudyDesign <-  readRDS('./Data/model_study_design.rds')

# Load species taxa groups and sensitivity scores
source('Data/taxa_tolerance.R')
source('Data/taxa_groups.R')

# load functions
source('Functions/predict_gllvm.R')


# Create new data

Variables = c("Mean", "Trend2", "Residual", "Seasonal", "l_SUD_N")
newdata <- list()
n_gradient = 50 # number of equaly spaced points within each gradient (variable)

for (i in Variables) {
  kk <- data.frame(cbind(rep(seq(-2.5,2.5, length.out = n_gradient),length(row.names(mAb$lvs))), 0, 0, 0, 0), rep(row.names(mAb$lvs), each = n_gradient))
  names(kk)[1] <- i
  kk[,1] <- kk[,1]*sd(mod.mat[,i])
  names(kk)[-1] <- c(Variables[!(Variables %in% i)], 'site')
  newdata[[length(newdata) + 1]] <- kk 
}

# simulate
simulations <- lapply(newdata, function(x) {
  LV <- data.frame(mAb$lvs[x[,'site'],])
  newX <- x[,1:5]
  
  if (names(x)[1] != 'Mean') {
  for (site in unique(StudyDesign$Site)) {
    newX[,'Mean'] <- mod.mat[StudyDesign$Site == site, 'Mean'][1]
    }
  }
  
  if (names(x)[1] != 'l_SUD_N') {
    for (site in unique(StudyDesign$Site)) {
      newX[,'l_SUD_N'] <- mean(mod.mat[StudyDesign$Site == site, 'l_SUD_N'])
    }
  }
  
  com <- pred.gllvm(mAb, Xnew = newX, nsim = 1, LV = LV, sim_dist = F)
  x[,'var'] <- names(x)[1]
  x[,'val'] <- x[,1]
  data.frame(cbind(x, com))
})

sim_data <- do.call(rbind, simulations)
out <- list() 

# EXPECTED COMMUNITY COMPOSITION
for (v in unique(sim_data$var)) {
  print(v)
  for (s in unique(sim_data$site)) {
    print(s)
    kk <- subset(sim_data, site == s & var == v)
    
    # Pre-allocate vectors for adjacent steps only
    n_steps <- nrow(kk) - 1
    beta_bray <- numeric(n_steps)
    beta_grad <- numeric(n_steps)
    beta_bal <- numeric(n_steps)
    
    # Only loop through adjacent pairs (i and i+1)
    for (i in 1:n_steps) {
      j <- i + 1
      kk_beta <- betapart::beta.pair.abund(log1p(kk[c(i,j), -c(1:9)]))
      
      beta_bray[i] <- unlist(kk_beta$beta.bray)
      beta_grad[i] <- unlist(kk_beta$beta.bray.gra)
      beta_bal[i] <- unlist(kk_beta$beta.bray.bal)
    }
    
    out_kk <-  data.frame(var = v, 
                          site = s,
                          val_i = kk$val[1:n_steps],
                          val_j = kk$val[2:(n_steps+1)],
                          beta_bray = beta_bray,
                          beta_grad = beta_grad,
                          beta_bal = beta_bal)
    out[[length(out) + 1]] <- out_kk
  }
}

out_data <- data.frame(do.call(rbind, out))

# RATE OF CHANGE
out_data$bray <- out_data$beta_bray / abs(out_data$val_i - out_data$val_j)
out_data$balanced <- out_data$beta_bal / abs(out_data$val_i - out_data$val_j)
out_data$gradient <- out_data$beta_grad / abs(out_data$val_i - out_data$val_j)

# PLOT RATE OF CHANGE
out_data %>% na.omit() %>% pivot_longer(c(bray, balanced, gradient), names_to = 'partition', values_to = 'rate') %>%
  mutate(
    partition = factor(partition, levels = c('bray', 'balanced', 'gradient'), labels = c('Bray-Curtis', 'Balanced', 'Gradient')),
    var = recode(var, Mean = ' mean', Trend2 = 'slow',  Seasonal = 'seasonal', Residual = 'fast', l_SUD_N = 'N eq. SUD')) %>%
  filter(var %in% c('slow', 'seasonal', 'fast')) %>%
  group_by(var, site, partition) %>%
  summarise(rate = mean(rate)) %>%
  ggplot(aes(x = partition, y = rate)) +
  geom_violin(aes(fill = var, group = interaction(partition, var)), bounds = c(0,1), colour = 'transparent', scale = 'width', alpha = 0.2, adjust = 1, trim = FALSE) +
  geom_jitter(height = 0, width = 0.2, shape = 16, aes(colour = var), size = 0.5, alpha = 0.5) +
  stat_summary(geom = 'linerange', aes(colour = var, group = interaction(partition, var)), position = position_dodge2(width = 0.9), fun.max = function(x) mean(x) + sd(x), fun.min = function(x) mean(x) - sd(x) ) +
  stat_summary(geom = 'point', shape = 21, fill = 'white',aes(colour = var, group = interaction(partition, var)), position = position_dodge2(width = 0.9), fun = function(x) mean(x)) +
  xlab('dissimilarity partition') +
  facet_wrap(~var, ncol = 3) + 
  theme(axis.ticks.x = element_blank(), 
        strip.background = element_blank(),
        strip.text  = element_text(face = 2),
        aspect.ratio = 1,
        legend.position = '') + 
  scale_colour_manual(values = c('grey20', 'coral', 'steelblue')) + 
  scale_fill_manual(values = c('grey20', 'coral', 'steelblue')) + 
  scale_y_continuous('rate of community change \n (Bray-Curtis / log-flow)', expand = c(0,0), limits = c(0,0.42))

ggsave('Figs/rate_com_change.png', dpi = 600, width = 18, height = 8, units = 'cm')
ggsave('Figs/rate_com_change.pdf', dpi = 600, width = 18, height = 8, units = 'cm')

