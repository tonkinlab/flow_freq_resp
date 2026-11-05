library(ggplot2)
library(cowplot)
library(tidyverse)
library(psych)

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

group_colours =   c('purple4','steelblue', 'darkolivegreen4', 'khaki3', 'gold2','darkorange', 'coral2', 'brown','grey20','grey20','grey20')

coef_data %>% 
  filter(var2 %in% c('fast', 'slow')) %>%
  mutate(
    group = factor(taxa_groups[sp], levels = c('Ephemeroptera', 'Plecoptera', 'Trichoptera', 'Diptera', 'Coleoptera', 'Gastropoda','Crustacea', 'Worms', 'Arachnida', 'Hemiptera', 'Megaloptera')),
    var2_label = recode(var2, 
                        fast = "atop(bold('Fast'), 'MAE: 0.70')",
                        slow = "atop(bold('Slow'), 'MAE: 1.68')")
  ) %>%
  ggplot(aes(x = coef, y = sp)) +
  facet_grid(group ~ var2_label, scales = 'free_y', space = 'free_y', labeller = label_parsed) + 
  theme(
    legend.position = 'bottom',
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.5),
    axis.line.x = element_line(linewidth = 0.5),
    legend.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(colour = 'grey95', linewidth = 0.2),
    panel.spacing = unit(0.5, "lines"),
    strip.background = element_blank(),
    strip.text.x = element_text(colour = 'grey20', size = 10, face = "plain", lineheight = 1.1),
    strip.text.y = element_text(angle = 0, face = "plain", size = 8, hjust = 0),
    panel.border = element_rect(colour = NA, fill = NA),
    panel.background = element_rect(fill = 'white', colour = NA)
  ) +
  geom_vline(xintercept = 0, lty = 2, col = 'grey50', linewidth = 0.2) +
  geom_segment(aes(xend = 0, x = coef, colour = group, alpha = sig), linewidth = 0.4) +
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
library(tidyverse)
library(ggrepel)

pred_x <- "Trend2"; pred_y <- "Residual"
level  <- 0.9
min_prev <- 0.1                                  # leave out rare taxa (huge ellipses)

V  <- vcov(mAb)                                   # covariance of all estimated parameters
ib <- which(rownames(V) == "b")                   # the species-specific coefficients are all labelled "b"
se_b <- sqrt(diag(V))[ib]
Xs <- mAb$sd$Xcoef                                # SEs as a taxa x predictors matrix
Bs <- mAb$params$Xcoef                            # estimates

# The "b" names don't say which entry is which taxon / predictor, so find each one by matching
# its standard error against mAb$sd$Xcoef (exact match up to rounding).
find_idx <- function(j, k) {
  d <- abs(se_b - Xs[j, k])
  if (min(d) > 1e-6 * max(Xs[j, k], 1e-12)) stop("no matching SE for ", rownames(Xs)[j], " / ", k)
  if (sum(d < 1e-9 * max(Xs[j, k], 1e-12)) > 1) warning("ambiguous SE match for ", rownames(Xs)[j], " / ", k)
  ib[which.min(d)]
}

ellipse_xy <- function(mu, S, level = 0.95, n = 100) {
  th <- seq(0, 2 * pi, length.out = n)
  r  <- sqrt(qchisq(level, df = 2))
  t(mu + r * t(chol(S)) %*% rbind(cos(th), sin(th)))   # n x 2
}

taxa <- rownames(Bs)
res <- map(seq_along(taxa), function(j) {
  i1 <- find_idx(j, pred_x); i2 <- find_idx(j, pred_y)
  S  <- V[c(i1, i2), c(i1, i2)]
  list(pt  = tibble(taxon = taxa[j], x = Bs[j, pred_x], y = Bs[j, pred_y],
                    r = S[1, 2] / sqrt(S[1, 1] * S[2, 2])),     # correlation between the two estimates
       ell = as_tibble(setNames(as.data.frame(ellipse_xy(c(Bs[j, pred_x], Bs[j, pred_y]), S, level)),
                                c("x", "y"))) %>% mutate(taxon = taxa[j]))
})
pts  <- bind_rows(map(res, "pt"))  %>% left_join(tibble(taxon = colnames(y), prev = colMeans(y > 0)), by = "taxon")
ells <- bind_rows(map(res, "ell")) %>% left_join(select(pts, taxon, prev), by = "taxon")

keep_t <- pts$taxon[pts$prev > min_prev]

# colour by taxonomic group (same groups and colours as the coefficient plot above)
group_levels <- c('Ephemeroptera', 'Plecoptera', 'Trichoptera', 'Diptera', 'Coleoptera', 'Gastropoda',
                  'Crustacea', 'Worms', 'Arachnida', 'Hemiptera', 'Megaloptera')
pts  <- pts  %>% mutate(group = factor(taxa_groups[taxon], levels = group_levels))
ells <- ells %>% mutate(group = factor(taxa_groups[taxon], levels = group_levels))
stopifnot(!anyNA(pts$group[pts$taxon %in% keep_t]))   # every plotted taxon needs an entry in taxa_groups

p <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_path(data = filter(ells, taxon %in% keep_t), aes(y, x, group = taxon, colour = group), linewidth = 0.3, , alpha = 0.1) +
  geom_point(data = filter(pts, taxon %in% keep_t), aes(y, x, colour = group), size = 0.8) +
  geom_text_repel(data = filter(pts, taxon %in% keep_t), aes(y, x, label = taxon, colour = group),
                  size = 2.3, max.overlaps = Inf, segment.size = 0.2) +
  scale_colour_manual("Group", values = setNames(group_colours, group_levels)) +
  scale_alpha_continuous("Prevalence", range = c(0.35, 1)) +      # rarer taxa are fainter
  labs(x = paste(pred_x, "coefficient"), y = paste(pred_y, "coefficient"),
       caption = paste0(100 * level, "% confidence ellipses, taxa with prevalence > ", min_prev)) +
  theme_bw() + theme(panel.grid = element_blank()) +
  coord_equal()
print(p)


