# LOAD PACKAGES  ----
library(gllvm)
library(tidyverse)
library(lubridate)
library(corrplot)
library(doParallel)
library(foreach)
library(ggplot2)
library(ggtext)
library(cowplot)

# Set cores to use for optimization
n_cores <- detectCores(logical = T) - 2

# Plotting preferences
theme_set(theme_bw())
theme_update(text = element_text(size = 8),
             panel.grid = element_blank())

# Set seed
set.seed(2)

outDir = file.path("./Models")
if (!dir.exists(outDir)) dir.create(outDir)

# PREP DATA ----
data = readRDS("./Data/Data_models.RDS")
Y = data$Ymod          # Species data
XData = data$EnvData   # Environmental data

# Transform env data
XData$Trend2 <- XData$Trend - XData$meanFlow
XData$l_SUD_N <- log1p(XData$SUD_N)
XData$SUD_N <- NULL
XData$Site <- as.factor(XData$Site)

#Trait data 
TrData = read.csv("Data/traits.csv")
rownames <- TrData$taxa_match

# Drop species not in abundance data  and select Life history + mobility traits
colnames(TrData)
TrData <- TrData[, c(2,10:21, 39:41, 42:45, 66:68)]

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

# create a new category for size (values based on Macro-invertebrate trait database)
TrData$SIZE <- (TrData$Size_1 * 2.5 + TrData$Size_2 * 7.5 + TrData$Size_3 * 15 + TrData$Size_4 * 30 + TrData$Size_5 * 50) / (TrData$Size_1 + TrData$Size_2 + TrData$Size_3 + TrData$Size_4 + TrData$Size_5)
TrData$l_SIZE <- log(TrData$SIZE)
TrData <- dplyr::select(TrData, l_SIZE, aq_ADUANDLAR2, aq_ADUORLAR2, DissPotential_LOW2, DissPotential_HIGH2, Mob_BURROWER2,Mob_SWIMMER2,Mob_ATTACHED2,CycPerYear_SEMI2, CycPerYear_PLURIV2)
 
names(TrData)
row.names(TrData) <- rownames

# Study design
StudyDesign <- data.frame(
  Year = as.numeric(year(as.Date(data$Spt$Date))),
  lat = data$Spt$lat,
  lon = data$Spt$lon,
  Site = (as.factor(data$Spt$Site)),
  Year_site = as.factor(paste0(data$Spt$Date,data$Spt$Site)))

#Create model matrix
mod.mat <- data.frame(model.matrix(~ meanFlow + Trend2 + Residual + Seasonal + l_SUD_N, data = XData)[,-1])
names(mod.mat) <- c('Mean' , 'Trend2', 'Residual','Seasonal', 'l_SUD_N')
mod.mat <- data.frame(scale(mod.mat,scale = FALSE)) # Maintain units so coefficients can be compared
mod.mat$int <- 1 # Add intercept to taxa-specific random effects

TrData <- na.omit(TrData)
y <- Y[, rownames(TrData)] 
y <- y[, order(colSums(y > 0), decreasing = T)]
TR <- TrData[names(y),]

saveRDS(TrData, 'Data/TrData_model.rds')
saveRDS(mod.mat, 'Data/model_matrix_env.rds')
saveRDS(StudyDesign, 'Data/model_study_design.rds')

# FIT MODELS  ----
cl <- makeCluster(n_cores)
registerDoParallel(cl)

out_ab <- foreach(i = 1:10, .packages = "gllvm") %dopar% {
  gllvm(y = y,
             X = data.frame(mod.mat),
             formula = ~ Mean + Trend2 +  Seasonal + Residual + l_SUD_N ,
             studyDesign = StudyDesign,
             control.start = list(n.init = 1, jitter.var = 0.01), # random start
             row.eff =  ~ corAR1(1|Year),
             lvCor = ~ (1|Site),
             control = list(TMB = TRUE,
                            trace = TRUE,
                            max.iter = 200000,
                            maxit = 200000,
                            optimizer = 'optim'), # options for optimiser
             trace = TRUE, # print progress
             num.lv = 2,
             family = 'negative.binomial',
             seed = primes::primes[i],
             method = 'VA')
}

stopCluster(cl)

### SELECT BEST MODEL  ###
logL <- lapply(out_ab, function(x) x[['logL']])
mAb <- out_ab[[which.max(logL)]]
saveRDS(mAb, paste0(outDir, '/mAb.rds'))

# With traits
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Fit models in parallel
out1 <- foreach(i = 1:10, .packages = "gllvm") %dopar% {
    gllvm(y = y,
           TR = TR,
           X = mod.mat,
           formula = ~ (Mean + Trend2 + Seasonal + Residual + l_SUD_N)*(l_SIZE + aq_ADUANDLAR2 + aq_ADUORLAR2 + DissPotential_LOW2 + DissPotential_HIGH2 + Mob_BURROWER2 + Mob_SWIMMER2 + Mob_ATTACHED2 + CycPerYear_SEMI2 + CycPerYear_PLURIV2),
           studyDesign = StudyDesign,
           lvCor = ~ (1|Site),
           num.lv = 2,
           seed = primes::primes[i],
           control.start = list(n.init = 1, jitter.var = 0.01), # random start
           beta0com = TRUE,
           control = list(TMB = TRUE,
                          trace = TRUE,
                          max.iter = 200000,
                          maxit = 200000,
                          optimizer = 'optim'), # options for optimizer
           randomX = ~ Mean + Trend2 + Seasonal + Residual + l_SUD_N + int,
           family = 'negative.binomial')
}

stopCluster(cl)
logL <- lapply(out1, function(x) x[['logL']])
mTr <- out1[[which.max(logL)]]
saveRDS(mTr, paste0(outDir, '/mTr.rds'))



