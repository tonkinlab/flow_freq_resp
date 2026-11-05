# LOAD PACKAGES  ----
library(gllvm)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(cowplot)

source('Functions/Interactive_flow_decomp.R')   # half_gaussian_kernel() and decompose_flow()


slow_vals   <- c(seq(200, 365*5, 100))      # slow kernel sizes (days), tried with the best single kernel as the fast kernel
single_vals <- c(30,40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 350, 400, 450, 500, 550, 600)   # refine around 90 so the one-kernel model gets its best shot
single_fast <- 1                            # dummy fast kernel for the one-kernel models
fast_max    <- 120                          # the best single kernel becomes the fast kernel; it must not exceed this
kernel_window <- max(slow_vals, single_vals)   # every slow / single kernel is evaluated over this many days, so a kernel of
# any size decays to ~0 instead of being cut off at its own size
seasonal_comp <- FALSE

# Plotting preferences
theme_set(theme_bw())
theme_update(text = element_text(size = 8),
             panel.grid = element_blank())

# Set seed
set.seed(2)

outDir = file.path("./Models")
if (!dir.exists(outDir)) dir.create(outDir)

# PREP DATA
data = readRDS("./Data/Data_models.RDS")
Y = data$Ymod          # Species data
# Tanytarsini is a tribe within the subfamily Chironominae. Merged for consistency within the dataset.
Y$Chironominae <- Y$Tanytarsini + Y$Chironominae
Y$Tanytarsini <- NULL
XData = data$EnvData   # Environmental data

# FLOW at the sampling dates ----
imp <- read.csv('Data/Flow/FlowData_2022_imputed.csv')   # imputed log flow, all sites and days
imp$Date <- as.Date(imp$Date)
obs <- data.frame(site = as.character(data$Spt$Site),    # one row per sample, same order as Y
                  Date = as.Date(data$Spt$Date))

# decomposition of the flow for one pair of kernels, at the sampling dates (same row order as Y)
decomp_obs <- function(slow, fast) {
  d <- map_dfr(unique(obs$site), function(s)
    suppressWarnings(decompose_flow(imp, c(slow = slow, fast = fast), seasonal_comp,
                                    s, unique(obs$Date[obs$site == s]),
                                    slow_window = kernel_window, fast_window = 2 * fast)))
  d[match(paste(obs$site, obs$Date), paste(d$site, d$Date)), ]
}

# Samples without enough flow history (burn-in) have no components. Shorter kernels need less history, so
# dropping the samples that are missing for the LONGEST kernels, once, keeps every combination on the same rows.
dec  <- decomp_obs(kernel_window, fast_max)
keep <- !is.na(dec$slow) & !is.na(dec$fast)
message(sum(!keep), ' of ', length(keep), ' samples dropped (not enough flow history for the longest kernels)')
Y <- Y[keep, ]
XData <- XData[keep, ]
data$Spt <- data$Spt[keep, ]
obs <- obs[keep, ]
dec <- dec[keep, ]
rownames(Y) <- NULL
rownames(XData) <- NULL
rownames(data$Spt) <- NULL
rownames(obs) <- NULL

XData$meanFlow <- dec$mean   # site mean log flow: does not depend on the kernels

# Transform env data
XData$l_SUD_N <- log1p(XData$SUD_N)
XData$SUD_N <- NULL
XData$Site <- as.factor(XData$Site)

#Trait data
TrData = read.csv("Data/traits.csv")
rownames(TrData) <- TrData$taxa_match

# Drop species not in abundance data and select Life history + mobility traits
colnames(TrData)
TrData <- TrData[, c(2,10:21, 39:41, 42:45, 66:68)]
TrData <- TrData[which(rownames(TrData) %in% names(Y)),]

# FUZZY-CODED TRAITS are modified so they add to one, and one category becomes the reference
# Larva and pupa (aq_LARANDPUP2) becomes reference
TrData$aq_ADUANDLAR2 <- TrData$aq_ADUANDLAR / (TrData$aq_ADUANDLAR + TrData$aq_ADUORLAR + TrData$aq_LARANDPUP)
TrData$aq_LARANDPUP2 <- TrData$aq_LARANDPUP / (TrData$aq_ADUANDLAR + TrData$aq_ADUORLAR + TrData$aq_LARANDPUP)
TrData$aq_ADUORLAR2 <- TrData$aq_ADUORLAR / (TrData$aq_ADUANDLAR + TrData$aq_ADUORLAR + TrData$aq_LARANDPUP)
TrData$aq_LARANDPUP2 <- NULL

# dissemination medium (DissPotential_MEDIUM2) becomes reference
TrData$DissPotential_LOW2 <-  TrData$DissPotential_LOW / (TrData$DissPotential_LOW + TrData$DissPotential_MEDIUM + TrData$DissPotential_HIGH)
TrData$DissPotential_MEDIUM2 <-  TrData$DissPotential_MEDIUM / (TrData$DissPotential_LOW + TrData$DissPotential_MEDIUM + TrData$DissPotential_HIGH)
TrData$DissPotential_HIGH2 <-  TrData$DissPotential_HIGH / (TrData$DissPotential_LOW + TrData$DissPotential_MEDIUM + TrData$DissPotential_HIGH)
TrData$DissPotential_MEDIUM2 <- NULL

# mobility (mob) crawler becomes reference
TrData$Mob_CRAWLER2 <- TrData$Mob_CRAWLER / (TrData$Mob_CRAWLER + TrData$Mob_BURROWER + TrData$Mob_SWIMMER + TrData$Mob_ATTACHED)
TrData$Mob_BURROWER2 <- TrData$Mob_BURROWER / (TrData$Mob_CRAWLER + TrData$Mob_BURROWER + TrData$Mob_SWIMMER + TrData$Mob_ATTACHED)
TrData$Mob_SWIMMER2 <- TrData$Mob_SWIMMER / (TrData$Mob_CRAWLER + TrData$Mob_BURROWER + TrData$Mob_SWIMMER + TrData$Mob_ATTACHED)
TrData$Mob_ATTACHED2 <- TrData$Mob_ATTACHED / (TrData$Mob_CRAWLER + TrData$Mob_BURROWER + TrData$Mob_SWIMMER + TrData$Mob_ATTACHED)
TrData$Mob_CRAWLER2 <- NULL

# maximum number of reproductive events, make UNIV reference
TrData$CycPerYear_UNIV2 <- TrData$CycPerYear_UNIV / (TrData$CycPerYear_UNIV + TrData$CycPerYear_SEMI + TrData$CycPerYear_PLURIV)
TrData$CycPerYear_SEMI2 <- TrData$CycPerYear_SEMI / (TrData$CycPerYear_UNIV + TrData$CycPerYear_SEMI + TrData$CycPerYear_PLURIV)
TrData$CycPerYear_PLURIV2 <- TrData$CycPerYear_PLURIV / (TrData$CycPerYear_UNIV + TrData$CycPerYear_SEMI + TrData$CycPerYear_PLURIV)
TrData$CycPerYear_UNIV2 <- NULL

# create a new variable for size (values based on Macro-invertebrate trait database)
TrData$SIZE <- (TrData$Size_1 * 2.5 + TrData$Size_2 * 7.5 + TrData$Size_3 * 15 + TrData$Size_4 * 30 + TrData$Size_5 * 50) / (TrData$Size_1 + TrData$Size_2 + TrData$Size_3 + TrData$Size_4 + TrData$Size_5)
TrData$l_SIZE <- log(TrData$SIZE)
TrData <- dplyr::select(TrData, l_SIZE, aq_ADUANDLAR2, aq_ADUORLAR2, DissPotential_LOW2, DissPotential_HIGH2, Mob_BURROWER2,Mob_SWIMMER2,Mob_ATTACHED2,CycPerYear_SEMI2, CycPerYear_PLURIV2)

names(TrData)

# Study design
StudyDesign <- data.frame(
  Year = as.numeric(year(as.Date(data$Spt$Date))),
  lat = data$Spt$lat,
  lon = data$Spt$lon,
  Site = (as.factor(data$Spt$Site)),
  Year_site = as.factor(paste0(data$Spt$Date,data$Spt$Site)))

mod.mat <- data.frame(Mean = XData$meanFlow,
                      l_SUD_N = XData$l_SUD_N)
stopifnot(!anyNA(mod.mat))

# flow components for one pair of kernels added to the model matrix (fast = NULL: one-kernel model, slow only)
add_flow <- function(slow, fast = NULL) {
  d  <- decomp_obs(slow, if (is.null(fast)) single_fast else fast)
  mm <- mod.mat
  mm$Trend2 <- d$slow - d$mean                    # slow component (deviation from the site mean)
  if (!is.null(fast)) mm$Residual <- d$fast       # fast component
  stopifnot(!anyNA(mm))
  mm <- data.frame(scale(mm, scale = FALSE))      # center but keep scale so coefficients can be compared
  mm$int <- 1                                     # Add intercept to taxa-specific random effects
  mm
}

# predictors used in the formulas below
preds_flow  <- c('Trend2', 'Residual')
preds_other <- c('Mean', 'l_SUD_N')
preds <- c(preds_flow, preds_other)

y <- Y[,-c(1,2)] # remove non-abundance columns
y_ab <- y[colSums(y>0)/nrow(y) > 0.1] # only occurrences > 0.1 for the no-traits model
y_ab <- y_ab[, order(colSums(y_ab > 0), decreasing = T)]




# FIT MODELS  ----
TMB::openmp(5, DLL = 'gllvm')

# abundance model, species-specific effects of every predictor (full rank) + unconstrained LVs with a random site effect
fit_ab <- function(mm, flow_terms, start = NULL) {
  gllvm(y = y_ab,
        X = data.frame(mm),
        formula = reformulate(c(flow_terms, preds_other)),
        num.lv = 2,
        studyDesign = StudyDesign,
        lvCor = ~ (1|Site),
        control.start = if (is.null(start)) list(n.init = 3) else list(start.fit = start),   # start from a previous fit
        control = list(TMB = TRUE, trace = TRUE, optimizer.trace = 1, maxit = 200000),
        family = 'negative.binomial')
}

# KERNEL SELECTION ----
res <- data.frame()                                   # one row per fitted model

# 1. single-kernel models (Trend2 only): all sizes, each starting from the previous fit
prev <- NULL
for (s in single_vals) {
  fit <- tryCatch(fit_ab(add_flow(s), 'Trend2', prev),
                  error = function(e) tryCatch(fit_ab(add_flow(s), 'Trend2'), error = function(e2) NULL))   # retry from scratch
  if (is.null(fit)) { message('failed: single ', s); next }
  prev <- fit
  res <- rbind(res, data.frame(type = 'one', slow = s, fast = NA_real_,
                               logL = fit$logL, df = attr(logLik(fit), 'df')))
  saveRDS(res, file.path(outDir, 'kernel_sweep_window.rds'))                        # progress is kept if R crashes
  message('single ', s, ': logL ', round(fit$logL, 1))
}
one <- res[res$type == 'one', ]
fast_best <- one$slow[which.max(one$logL)]                                   # best single kernel
message('best single kernel: ', fast_best)
if (fast_best > fast_max) stop('best single kernel (', fast_best, ') is longer than fast_max (', fast_max, '): raise fast_max and rerun')

prev <- NULL
for (s in slow_vals) {
  print(s)
  if (!any(res$type == 'two' & res$slow == s)) {                              # skip sizes already fitted (the one-kernel rows do not count)
    fit <- tryCatch(fit_ab(add_flow(s, fast_best), preds_flow, prev),   # start from the previous slow fit
                    error = function(e) tryCatch(fit_ab(add_flow(s, fast_best), preds_flow), error = function(e2) NULL))
    if (is.null(fit)) { message('failed: slow ', s, ' fast ', fast_best); next }
    prev <- fit
    res <- rbind(res, data.frame(type = 'two', slow = s, fast = fast_best,
                                 logL = fit$logL, df = attr(logLik(fit), 'df')))
    saveRDS(res, file.path(outDir, 'kernel_sweep_window.rds'))
    message('slow ', s, ' / fast ', fast_best, ': logL ', round(fit$logL, 1))
  }
}

sweep <- res                   
sweep <- bind_rows(sweep) %>%
  mutate(AIC = -2 * logL + 2 * df, BIC = -2 * logL + log(nrow(y_ab)) * df,
         dAIC = AIC - min(AIC), dBIC = BIC - min(BIC))

aic_one  <- min(sweep$AIC[sweep$type == 'one'])                            # best single-kernel model
size_one <- sweep$slow[sweep$type == 'one'][which.min(sweep$AIC[sweep$type == 'one'])]

one <- sweep %>% filter(type == 'one')
two <- sweep %>% filter(type == 'two')
col_one <- scale_colour_viridis_c('AIC', direction = -1, limits = range(one$AIC))   # best (lowest AIC) = brightest
col_two <- scale_colour_viridis_c('AIC', direction = -1, limits = range(two$AIC))

# top row: AIC against kernel size
p_one <- ggplot(one, aes(-slow, AIC)) +
  geom_line(colour = 'grey60') + geom_point(aes(colour = AIC), size = 2) + col_one +
  labs(x = 'Single kernel (days)', y = 'AIC') +
  theme(legend.position = 'none')

p_two <- ggplot(two, aes(-slow/365, AIC)) +
  geom_line(colour = 'grey60') + geom_point(aes(colour = AIC), size = 2) + col_two +
  labs(x = 'Slow kernel (years)', y = 'AIC') +
  theme(legend.position = 'none')

# bottom row: the kernels themselves (weight of each day back from the sampling date, 1 = weight of the day itself),
# coloured by the AIC of the model that used them
kern_one <- map2_dfr(one$slow, one$AIC, function(s, a) {
  k <- half_gaussian_kernel(s / 2, kernel_window); lag <- 0:kernel_window
  tibble(size = s, AIC = a, lag = lag, w = k / k[1]) %>% filter(lag <= 2 * s) })
kern_two <- map2_dfr(two$slow, two$AIC, function(s, a) {
  k <- half_gaussian_kernel(s / 2, kernel_window); lag <- 0:kernel_window
  tibble(size = s, AIC = a, lag = lag, w = k / k[1]) %>% filter(lag <= 2 * s) })

p_kern_one <- ggplot(kern_one, aes(-lag, w, group = size, colour = AIC)) +
  geom_line(linewidth = 0.6) + col_one +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(-500, 0)) +  
  labs(x = 'Days before the sampling date', y = 'Relative kernel weight')

p_kern_two <- ggplot(kern_two, aes(-lag/365, w, group = size, colour = AIC)) +
  geom_line(linewidth = 0.6) + col_two +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(-5, 0)) +
  labs(x = 'years before the sampling date', y = 'Relative kernel weight')

plot_grid(p_one, p_two, p_kern_one, p_kern_two, nrow = 2, align = 'v', axis = 'lr', rel_heights = c(2, 3))



print(sweep %>% arrange(AIC) %>% head(10) %>% as.data.frame(), digits = 6)
best2 <- sweep %>% filter(type == 'two') %>% slice_min(AIC, n = 1, with_ties = FALSE)
best1 <- sweep %>% filter(type == 'one') %>% slice_min(AIC, n = 1, with_ties = FALSE)
cat(sprintf('Best two kernels: slow %d / fast %d (AIC %.1f). Best single kernel: %d (AIC %.1f). dAIC = %.1f, dBIC = %.1f\n',
            best2$slow, best2$fast, best2$AIC, best1$slow, best1$AIC, best1$AIC - best2$AIC, best1$BIC - best2$BIC))

# Does each component add something? Likelihood-ratio tests inside the best two-kernel model (nested models)

# Save intermediate data for use in other scripts
saveRDS(TrData, 'Data/TrData_model.rds')
saveRDS(mod.mat, 'Data/model_matrix_env.rds')
saveRDS(StudyDesign, 'Data/model_study_design.rds')

saveRDS(mAb, paste0(outDir, '/mAb.rds'))


mod.mat <- add_flow(900, 80)

mAb <- fit_ab(mod.mat, preds_flow)


# Fit model
TrData <- na.omit(TrData)
y_traits <- y[, rownames(TrData)]
y_traits <- y_traits[, order(colSums(y_traits > 0), decreasing = T)]
TR <- TrData[names(y_traits),]

mTr <- gllvm(y = y_traits,
             TR = TR,
             X = mod.mat,
             formula = as.formula(paste0('~ (', paste(preds, collapse = ' + '), ') * (',
                                         paste(names(TrData), collapse = ' + '), ')')),
             studyDesign = StudyDesign,
             lvCor = ~ (1|Site),
             num.lv = 2,
             beta0com = FALSE,
             control = list(TMB = TRUE, trace = TRUE, optimizer.trace = 1, maxit = 200000),
             randomX = reformulate(preds),  
             family = 'negative.binomial')

saveRDS(mTr, paste0(outDir, '/mTr.rds'))

