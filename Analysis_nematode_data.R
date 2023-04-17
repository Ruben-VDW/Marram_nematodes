# Script to analyse nematode functional group data from marram grass
# Made by Ruben Van De Walle

# Part one: selecting environmental variables 
#  best explaining the variation within the data
################################################

## Load packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(Hmsc)
library(mcmcplots)
library(ggforce)
library(ggExtra)
library(corrplot)

#setwd("~")
setwd("~/ENDURE/Results/Biodiveristy/Soil community/uploaden")

## importing data
nematodes <- read.csv("marram_grass_nematodes_data.csv")
nematodes <- subset(nematodes, select = c(-X.1))
nematodes$country <- factor(nematodes$country)
nematodes$transect <- factor(nematodes$transect)
nematodes$district <- factor(nematodes$district)

################################################

## 1. Data visualization

# nematodes found in roots
Pgraphs <- nematodes
Pgraphs <- subset(Pgraphs, select=c(country, site, Bact_R, Fung_R, Plant_R, Omn_R, Pred_R, Unk_R)) %>%
  pivot_longer(., cols = c(Bact_R, Fung_R, Plant_R, Omn_R, Pred_R, Unk_R), names_to = "niche", values_to = "Count")

p <- ggplot(Pgraphs, aes(x=niche, y=Count, fill = niche)) + 
  geom_boxplot() + theme_classic() +
  labs(x="", y = "Nematodes per g dry roots") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        legend.position="none") #+ ylim(c(0,1000))
p + geom_hline(yintercept=1000, linetype="dashed")


# nematodes found in soil
Pgraphs <- subset(nematodes_no_NA, select=c(country, site, Bact_S, Fung_S, Plant_S, Omn_S, Pred_S, Unk_S)) %>%
  pivot_longer(., cols = c(Bact_S, Fung_S, Plant_S, Omn_S, Pred_S, Unk_S), names_to = "niche", values_to = "Count")


p2 <- ggplot(Pgraphs, aes(x=niche, y=Count, fill = niche)) + 
  geom_boxplot() + theme_classic() + 
  labs(x="", y = "Nematodes per g dry soil") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        legend.position="none") #+ ylim(c(0,10))
p2 + geom_hline(yintercept=10, linetype="dashed")


# Nematodes in roots per biogeographic region
Pgraphs <- nematodes
p3 <- ggplot(Pgraphs, aes(x=district, y=Total_R, fill = district)) +
  geom_boxplot() + theme_classic() +
  labs(x="", y = "Nematodes per g dry roots") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        legend.position = "none")
p3 


# Nematodes in roots per biogeographic region
p4 <- ggplot(Pgraphs, aes(x=district, y=Total_S, fill = district)) +
  geom_boxplot() + theme_classic() + 
  labs(x="", y = "Nematodes per g dry soil") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        legend.position="none") #+ ylim(c(0,30))
p4 + geom_hline(yintercept=30, linetype="dashed")

#####################################################



## 2. JSDM
## 2.1 Fitting the model

# Selecting data sets
#####################
nematodes_no_NA <- nematodes %>% drop_na(M10, Vitality) # dropping NA's
levels(nematodes_no_NA$transect)[12:13] <- c("ODK", "ODK")

summary(nematodes_no_NA$transect)
summary(nematodes_no_NA$country)
ggplot(nematodes_no_NA) + geom_bar(aes(x=transect))

# selecting explaining variables
env <- subset(nematodes_no_NA, select = c("P10", "M10", "Vitality", "senecio", "district", "transect", "sample_ID"))
rownames(env) <- nematodes_no_NA$sample_ID
XData <- env #this way, env still contains the raw data

# formatting species data
spp <- subset(nematodes_no_NA, select=c(Bact_S, Fung_S, Plant_S, Omn_S, Pred_S, Unk_S, 
                                        Bact_R, Fung_R, Plant_R,  Omn_R, Pred_R, Unk_R, sample_ID)) 
rownames(spp) <- spp$sample_ID
YData <- spp

# the data sets
head(XData)
head(YData)

# Make sure the rows of XData and YData match
XData <- XData %>% arrange(sample_ID)
YData <- YData %>% arrange(sample_ID)
Y <- as.matrix(YData[-13])
Y <- log(Y+1)


## Data preparation
###################
summary(XData)
XData$district <- factor(XData$district)
XData$senecio <- factor(XData$senecio)
XData$sample_ID <- factor(XData$sample_ID)

# We center all numerical environmental variables around their mean/median
XData$P10 <- XData$P10 - mean(XData$P10)
XData$M10 <- XData$M10 - mean(XData$M10) 
XData$Vitality <- XData$Vitality - median(XData$Vitality) # because this is a categorical variable that will be modelled as a continuous variable

# Set reference level to most common factor
XData$district <- relevel(XData$district, ref = "Flemish dunes")


## Model structure
###################
# We define the model
# Based on earlier data exploration, we assume the normal distribution as appropriate for log-transformed abundance data

# We consider all environmental covariates, and assume a quadratic effect for proportion and spatial configuration
XFormula_1 <- ~ P10 + M10 + Vitality + district + senecio

XFormula_2  <- ~ poly(P10, degree = 2, raw = TRUE) + M10 + Vitality + district + senecio

XFormula_3 <- ~ P10 + M10 + poly(Vitality, degree = 2, raw = TRUE) + district + senecio

XFormula_4 <- ~ poly(P10, degree = 2, raw = TRUE) + M10 + poly(Vitality, degree = 2, raw = TRUE) + district + senecio


XFormulas <- c(XFormula_1, XFormula_2, XFormula_3, XFormula_4)

# Implementing random effects
studyDesign_endure <- data.frame(sample = as.factor(XData$sample_ID), cluster = as.factor(XData$transect))
rL_soil_sample <- HmscRandomLevel(units = studyDesign_endure$sample)
rL_soil_transect <- HmscRandomLevel(units = studyDesign_endure$cluster)


## Model fitting
#################

# Set parameters
model.directory <- "models/"
nChains <- 4
thins <- c(100) #1000
niter <- c(1000)
i <- 1

# Run MCMC
for (thin in thins) {
  for (samples in niter) {
    for (XFormula in XFormulas){
      soil.setup <- Hmsc(Y=Y, XData=XData, XFormula=XFormula, studyDesign = studyDesign_endure,
                         ranLevels = list(sample = rL_soil_sample, cluster = rL_soil_transect), distr="normal")
      
      transient <- (samples/2)*thin # the number of MCMC steps that are executed before starting recording posterior samples
      verbose <- transient
      soil.jsdm <- sampleMcmc(soil.setup, nChains=nChains, thin=thin,
                              samples=samples, transient=transient,
                              verbose=verbose, nParallel=nChains)
      
      if (i == 1){rest <- "_PVit"}
      else if (i == 2){rest <- "_P2"}
      else if (i == 3){rest <- "_Vit2"}
      else if (i == 4){rest <- "_PVit2"}
      else {rest <- as.character(i)}
      filename <- paste0(model.directory, "small_soil_community_log_chains_", as.character(nChains),
                         "_samples_", as.character(samples),
                         "_thin_", as.character(thin), rest, "_a12_500")
      save(soil.jsdm, file=filename)
      i <- i + 1
    }
  }
}

######################################################



## 2.2 Evaluate convergence
load("models/soil_community_log_chains_4_samples_1000_thin_100_PVit")

## MCMC convergence
###################
# We first extract the posterior distribution from the model object and convert
# it into a coda object.
# Note that Hmsc uses coda only to examine convergence, whereas other operations
# (e.g. exploring of parameters or making predictions) are conducted straight
# from the Hmsc object soil.jsdm rather than the coda object.

soiljsdm_post <- convertToCodaObject(soil.jsdm)
#soiljsdm_post <- convertToCodaObject(root.jsdm)

## Examine convergencence in estimated parameters
# Which parameters were estimated?
names(soiljsdm_post)

# Visual assessment
mcmcplot(soiljsdm_post$Beta) # mcmcplots package
mcmcplot(soiljsdm_post$Gamma)
mcmcplot(soiljsdm_post$V)
mcmcplot(soiljsdm_post$Eta[[1]])
mcmcplot(soiljsdm_post$Psi[[1]])
mcmcplot(soiljsdm_post$Delta[[1]])
mcmcplot(soiljsdm_post$Lambda[[1]])
mcmcplot(soiljsdm_post$Omega[[1]])

# Formal diagnostics
# We first examine the effective size of the posterior sample.
# As there are three chains with a sample of 1000 from each, the actual sample 
# size is 3000. Thus, ideally the effective sample size would be 3000 as well.
# But in presence of autocorrelation, the effective sample size will be lower 
# than the actual sample size
ess.beta <- effectiveSize(soiljsdm_post$Beta)
hist(ess.beta, xlab = expression("Effective sample size" ~ beta ~ ""))

# We then examine the Gelman diagnostics, i.e. the Potential scale reduction factors
# This diagnostic compares if different chains (here we have 3 chains) give consistent results
# Ideally the value of this diagnostic would be close to one.
# As you increase thinning, you should observe the values getting closer to one.
# The Gelman diagnostic is often more informative diagnostic than the effective sample size
psrf.beta <- gelman.diag(soiljsdm_post$Beta, multivariate=FALSE)$psrf
summary(psrf.beta[, "Point est."])
summary(psrf.beta[, "Upper C.I."])

psrf.gamma <- gelman.diag(soiljsdm_post$Gamma, multivariate=FALSE)$psrf
summary(psrf.gamma[, "Point est."])
summary(psrf.gamma[, "Upper C.I."])

psrf.V <- gelman.diag(soiljsdm_post$V, multivariate=FALSE)$psrf
summary(psrf.V[, "Point est."])
summary(psrf.V[, "Upper C.I."])

psrf.Delta <- gelman.diag(soiljsdm_post$Delta[[1]], multivariate=FALSE)$psrf
summary(psrf.Delta[, "Point est."])
summary(psrf.Delta[, "Upper C.I."])

psrf.Psi <- gelman.diag(soiljsdm_post$Psi[[1]], multivariate=FALSE)$psrf
summary(psrf.Psi[, "Point est."])
summary(psrf.Psi[, "Upper C.I."])

psrf.Eta <- gelman.diag(soiljsdm_post$Eta[[1]], multivariate=FALSE)$psrf
summary(psrf.Eta[, "Point est."])
summary(psrf.Eta[, "Upper C.I."])

psrf.Lambda <- gelman.diag(soiljsdm_post$Lambda[[1]], multivariate=FALSE)$psrf
summary(psrf.Lambda[, "Point est."])
summary(psrf.Lambda[, "Upper C.I."])

psrf.Omega <- gelman.diag(soiljsdm_post$Omega[[1]], multivariate=FALSE)$psrf
summary(psrf.Omega[, "Point est."])
summary(psrf.Omega[, "Upper C.I."])

#making the plot
df_GR <- data.frame(GR_est = c(psrf.beta[, "Point est."], 
                               psrf.gamma[, "Point est."],
                               psrf.V[, "Point est."],
                               psrf.Eta[, "Point est."], 
                               psrf.Delta[, "Point est."],
                               psrf.Psi[, "Point est."],
                               psrf.Lambda[, "Point est."], 
                               psrf.Omega[, "Point est."]),
                    GR_CL = c(psrf.beta[, "Upper C.I."], 
                              psrf.gamma[, "Upper C.I."],
                              psrf.V[, "Upper C.I."],
                              psrf.Eta[, "Upper C.I."], 
                              psrf.Delta[, "Upper C.I."],
                              psrf.Psi[, "Upper C.I."],
                              psrf.Lambda[, "Upper C.I."], 
                              psrf.Omega[, "Upper C.I."]),
                    parameter = rep(c("beta", "mu", "V", "eta", "delta", "psi", 
                                      "lambda", "omega"), 
                                    c(length(psrf.beta[, 1]), 
                                      length(psrf.gamma[, 1]),
                                      length(psrf.V[, 1]),
                                      length(psrf.Eta[, 1]), 
                                      length(psrf.Delta[, 1]),
                                      length(psrf.Psi[, 1]),
                                      length(psrf.Lambda[, 1]),
                                      length(psrf.Omega[, 1]))))
df_GR$parameter <-factor(df_GR$parameter, levels = c('beta','mu', "V", "eta", 
                                                     "psi", "delta", "lambda", 
                                                     "omega"), ordered = TRUE)

cols <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "#CD9600", "#00BE67", 
          "#00A9FF", "#FF61CC")

# Point estimate
ggplot(df_GR, aes(x = parameter, y = GR_est)) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "firebrick") +
  geom_sina(aes(colour = parameter)) +
  labs(x = "", y = "Potential scale reduction factor") +
  #scale_y_continuous(breaks = seq(1, 1.006, by = 0.001)) +
  scale_x_discrete(labels = c(expression(beta), expression(mu), 'V', 
                              expression(eta), expression(psi), expression(delta),
                              expression(lambda), expression(omega))) +
  theme(axis.text.y =element_text(size=14),
        axis.text.x =element_text(size=22, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=16, face="bold"),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        #panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "") +
  scale_colour_manual(values = cols)

# Upper confidence limit
ggplot(df_GR, aes(x = parameter, y = GR_CL)) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "firebrick") +
  geom_sina(aes(colour = parameter)) +
  labs(x = "", y = "Potential scale reduction factor") +
  #scale_y_continuous(breaks = seq(1, 1.006, by = 0.001)) +
  scale_x_discrete(labels = c(expression(beta), expression(mu), 'V', 
                              expression(eta), expression(psi), expression(delta),
                              expression(lambda), expression(omega))) +
  theme(axis.text.y =element_text(size=14),
        axis.text.x =element_text(size=22, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title = element_text(size=16, face="bold"),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        #panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        legend.position = "") +
  scale_colour_manual(values = cols)

## Model fit
############
# We evaluate model fit here in terms of how well the model is able to discriminate abundances.
# This can be done by computing the model's explanatory power in terms of e.g. TjurR2 or AUC
# Below we make a plot that compares these two to each other.
# Note that the AUC values are much higher than the TjurR2 values, even if they measure 
# the fit of exactly the same model.
# Thus what is high or low value of model fit depends very much of the units used to measure it.

# RMSE measures how close the best predictions are to the data, i.e. a measure of accuracy
# R² values can be expected to be symmetrically distributed around zero if predictions are worse than by chance

# Explanatory power
pred <- computePredictedValues(soil.jsdm, nParallel = nChains)
MF_randENV <- evaluateModelFit(hM=soil.jsdm, predY=pred)
plot(MF_randENV$RMSE, MF_randENV$R2, main = "Explanatory power", xlab = "RMSE", ylab = "R2")

summary(MF_randENV$R2)
summary(MF_randENV$RMSE)

# Predictive power (2-fold cross-validation)
partition <- createPartition(soil.jsdm, nfolds = 2)
#predCV <- computePredictedValues(soil.jsdm, partition = partition, nParallel = nChains)
#save(predCV, file="predrandENV1_predictivePower")
# load("predrandENV1_predictivePower")
#MFCV_randENV <- evaluateModelFit(hM = soil.jsdm, predY = predCV)
#save(MFCV_randENV, file="modelfitrandENV1_predictivePower")
load("modelfitrandENV1_predictivePower")
plot(MFCV_randENV$RMSE, MFCV_randENV$R2, main = "Predictive power", xlab = "RMSE", ylab = "R2")


# The usual unconditional cross-validation can be used to assess
# the extent of overfitting of the fixed effects, which happens when the
# predictive power of the model is much lower than the explanatory power.
plot(MF_randENV$RMSE~MFCV_randENV$RMSE, xlab = "Predictive power", ylab = "Explanatory power", type = "n")
abline(a = 0, b = 1, col = "firebrick")
points(MF_randENV$RMSE~MFCV_randENV$RMSE, pch = 16, cex = 0.8)

df_RMSE <- data.frame("Explanatory power" = MF_randENV$RMSE, 
                      "Predictive power" = MFCV_randENV$RMSE)
p <- ggplot(df_RMSE, aes(x = Predictive.power, y = Explanatory.power)) +
  geom_hline(yintercept = 0, colour = "lightgrey", size = 1) + 
  geom_vline(xintercept = 0, colour = "lightgrey", size = 1) +
  geom_abline(intercept = 0, slope = 1, colour = "firebrick") +
  geom_point(size = 2) +
  scale_x_continuous("Predictive power (RMSE)", breaks = seq(0, 1.4, by = 0.1), 
                     limits = c(-0.01, 1.4)) +
  scale_y_continuous("Explanatory power (RMSE)", breaks = seq(0, 1.2, by = 0.1), 
                     limits = c(-0.01, 1.2)) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
ggMarginal(p, type = "boxplot", size = 12, 
           xparams = list(outlier.colour = "red"),
           yparams = list(outlier.colour = "red"))


plot(MF_randENV$R2~MFCV_randENV$R2, main = "R2", xlab = "Predictive power",
     ylab = "Explanatory power", type = "n")
abline(a = 0, b = 1, col = "firebrick")
points(MF_randENV$R2~MFCV_randENV$R2, pch = 16, cex = 0.8)

df_R2 <- data.frame("Explanatory power" = MF_randENV$R2, 
                    "Predictive power" = MFCV_randENV$R2)
p2 <- ggplot(df_R2, aes(x = Predictive.power, y = Explanatory.power)) +
  geom_abline(intercept = 0, slope = 1, colour = "firebrick") +
  geom_point(size = 2) +
  scale_x_continuous("Predictive power (R2)", breaks = seq(0, 0.3, by = 0.05), 
                     limits = c(0, 0.3)) +
  scale_y_continuous("Explanatory power (R2)", breaks = seq(0, 1, by = 0.05), 
                     limits = c(0, 1)) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))
ggMarginal(p2, type = "boxplot", size = 12, 
           xparams = list(outlier.colour = "red"),
           yparams = list(outlier.colour = "red"))

# Summary measures
summary(MF_randENV$RMSE)
summary(MFCV_randENV$RMSE)
summary(MF_randENV$R2)
summary(MFCV_randENV$R2)

#compute WAIC
computeWAIC(soil.jsdm)

######################################################



## 3.2 Results
# Prevalence
prev <- data.frame(P = colMeans(Y), species = colnames(Y))
prev <- prev %>% arrange(desc(P))
prev <- cbind(prev, index = 1:nrow(prev))


# Means and medians per covariate
means <- aggregate(subset(XData, select = c("P10", "M10")),
                   by = list(XData$district), mean)
median <- aggregate(subset(XData, select = c("Vitality")),
                    by = list(XData$district), median)

# We first examine the beta-parameters
# The heatmap shows whether the estimates are positive or negative with at least 0.95
# posterior probability.
postBeta <- getPostEstimate(soil.jsdm, parName="Beta")
plotBeta(soil.jsdm, post = postBeta, param = "Sign", plotTree = FALSE,
         supportLevel = 0.95, spNamesNumbers = c(TRUE, FALSE), covNamesNumbers = c(TRUE, FALSE))

# Make own plot
# Dataframe
betas <- c("intercept", "Prop", "Moran's I", "Vitality", "DIST_Boulonnais", "DIST_Renodunal", "DIST_UK", "DIST_Wadden", "Senecio")
support <- unlist(as.data.frame(postBeta$support), use.names = TRUE)
supportNeg <- unlist(as.data.frame(postBeta$supportNeg), use.names = TRUE)
support_df <- data.frame(supportPos = support, supportNeg = supportNeg, 
                         species = rep(colnames(soil.jsdm$Y), each = length(betas)))
values <- apply(support_df, 1, function(i) {
  if (i["supportNeg"] > 0.95) {
    return("-")
  } else if (i["supportPos"] > 0.95) {
    return("+")
  } else {
    return("0")
  }
})

beta_df <- data.frame(species = factor(rep(colnames(soil.jsdm$Y), each = length(betas)), 
                                       #levels = rev(prev$species),
                                       levels = c("Unk_S", "Pred_S", "Omn_S", "Fung_S", "Plant_S", "Bact_S",  
                                                  "Unk_R", "Pred_R", "Omn_R", "Fung_R", "Bact_R", "Plant_R"),
                                       ordered = TRUE),
                      beta = factor(rep(betas, length(colnames(soil.jsdm$Y))), 
                                    levels = betas, 
                                    ordered = TRUE),
                      Legend = factor(values, 
                                      levels = c("+", "0", "-"), 
                                      ordered = TRUE))
#levels(beta_df$species)

# Plot
# list(c(0, "#440154FF"), c(0.5, "73D055FF"), c(1, "#FDE725FF"))
beta_df <- subset(beta_df, beta != "intercept")
p <- ggplot(beta_df, aes(x=beta, y=species, fill=Legend)) + theme_bw() +
  geom_tile(colour = "black") +
  scale_fill_manual(values = c("#440154FF", "ghostwhite",  "#FDE725FF")) +
  labs(x = "", y = "")+
  theme(axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1,
                                   margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16, face="bold"))
p



## Variance partitioning
# How much of the variation in species abundance can be explained by a variable?

# Compute the variance partitioning
# Large vs small scale
VP_large <- computeVariancePartitioning(soil.jsdm,
                                        group=c(2, 1, 1, 1, 2, 2, 2, 2, 1), 
                                        groupnames = c("small scale", "large scale"))

soil.jsdm$covNames
# variance proportion for each group and species
# Mean over all species
apply(VP_large$vals, 1, summary)
round(apply(VP_large$vals, 1, summary), digits = 3)

# Fixed effects
VP <- computeVariancePartitioning(soil.jsdm, 
                                  group=c(1, 2, 3, 4, 1, 1, 1, 1, 5), 
                                  groupnames = c("District", "Proportion", "Configuration", "Vitality", "Senecio"))

# variance proportion for each group and species
# Mean over all species
apply(VP$vals, 1, summary)
round(apply(VP$vals, 1, summary), digits = 3)

boxplot(t(VP$vals), ylab = "Proportion of variance")

# Plant parasites in roots
VP$vals[, "Plant_R"]
round(VP$vals[, "Plant_R"], digits = 3)

# Plot with package
plotVariancePartitioning(soil.jsdm, VP = VP, 
                         cols = c("red", "blue", "green", "yellow", "purple", "grey",
                                  "pink", "orange"))

## Make own plot
# Dataframe
VP_temp <- as.data.frame(VP$vals)
species <- names(VP_temp)
components <- row.names(VP_temp)
VP_df <- data.frame(variance = unlist(VP_temp, use.names = TRUE), 
                    species = rep(species, each = length(components)),
                    Legend = rep(components, length(species)),
                    Scale = rep(c("Large scale", "Small scale", "Small scale", "Small scale", "Small scale",
                                  "Small scale", "Large scale"), length(species)))
VP_df$Legend <- factor(VP_df$Legend, levels = c("District", "Proportion", "Configuration", "Vitality", "Senecio",
                                                "Random: sample", "Random: cluster"), 
                       ordered = TRUE)
VP_df$species <- factor(VP_df$species, levels = c("Plant_R", "Bact_R", "Fung_R", "Omn_R", "Pred_R", "Unk_R",
                                                  "Plant_S", "Bact_S", "Fung_S", "Omn_S", "Pred_S", "Unk_S"), ordered = TRUE)
VP_df$Scale <- factor(VP_df$Scale, levels = c("Small scale", "Large scale"), 
                      ordered = TRUE)

# Colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n <- length(components)
cols <- gg_color_hue(n)
#cols <- c("coral", "forestgreen", "cyan3", "darkorange3", "grey", "yellow", "brown1")

# Plot
p <- ggplot(data = VP_df, aes(x=species, y=variance, fill=Legend, pattern = Scale)) +
  geom_bar_pattern(stat="identity", width = 1, colour = "black",
                   pattern_fill = "black",
                   pattern_density = 0.01,
                   pattern_spacing = 0.01,
                   pattern_key_scale_factor = 0.6) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand = c(0.01, 0)) +
  #scale_x_discrete(labels = c(1, rep("", 8), "10", rep("", 9), "20", rep("", 9),
  #                            "30", rep("", 9), "40", rep("", 9), 50)) +
  scale_fill_manual(values = cols) +
  scale_pattern_manual(values = c("Large scale" = "circle", "Small scale" = "none")) +
  labs(x = "", y = "Variance proportion", pattern = "") + 
  guides(fill = guide_legend(order = 1, override.aes = list(pattern = "none")),
         pattern = guide_legend(order = 2, override.aes = list(fill = "white"))) +
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, margin = margin(t = 5, r = 0, b = 0, l = 0), angle=90),
        axis.line = element_line(color = "gray40", size = 1, linetype = "solid"),
        axis.ticks.length = unit(5, "pt"),
        axis.title.y = element_text(size=16, face="bold", 
                                    margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size=16, face="bold", 
                                    margin = margin(t = 15, r = 0, b = 0, l = 0)),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16, face="bold"))
p

levels(VP_df$Legend)
levels(VP_df$Scale)
summary(subset(VP_df, Scale == "Small scale")$variance)
summary(subset(VP_df, Scale == "Large scale")$variance)


## Inference of biotic interactions?
# We next illustrate the species associations revealed by the random effects 
# with the corrplot function.
OmegaCor <- computeAssociations(soil.jsdm)
supportLevel <- 0.95
toPlot <- ((OmegaCor[[1]]$support > supportLevel) +
             (OmegaCor[[1]]$support < (1-supportLevel))>0) * OmegaCor[[1]]$mean
plotOrder <- corrMatOrder(OmegaCor[[1]]$mean, order="AOE")

# list(c(0, "#440154FF"), c(0.5, "73D055FF"), c(1, "#FDE725FF"))
corrplot(toPlot[plotOrder, plotOrder], method = "color", tl.cex = 0.8, tl.col = "black",
         col = colorRampPalette(c("#FDE725FF", "ghostwhite", "#440154FF"))(255),
         cl.cex = 1, cl.length = 11,
         addgrid.col = "lavender")

