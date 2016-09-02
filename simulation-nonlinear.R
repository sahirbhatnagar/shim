##################################
# R source code file for conducting NONLINEAR simulations
# on Mamouth cluster
# Git: this is on the eclust repo, sim2-modules-mammouth branch
# Created by Sahir,  August 31, 2016
# Updated: 
# Notes:
# note that I scp the files locally from the sim2-modules-mammouth branch to mammouth using:
# scp functions.R packages.R simulation.R simulation2.sh bhatnaga@cgreenwo-mp2.ccs.usherbrooke.ca:/home/bhatnaga/coexpression/august2016simulation/linear
# qsub command run from /home/bhatnaga/coexpression/august2016simulation/linear on mammouth:
# for i in 1 ; do qsub -v index=$i simulation2.sh ; done
##################################


rm(list=ls())
source("packages.R")
source("functions.R")

options(digits = 4, scipen = 999)

# source("/home/bhatnaga/coexpression/may2016simulation/sim2-modules-mammouth/packages.R")
# source("/home/bhatnaga/coexpression/may2016simulation/sim2-modules-mammouth/functions.R")
# source("/home/bhatnaga/coexpression/august2016simulation/linear/packages.R")
# source("/home/bhatnaga/coexpression/august2016simulation/linear/functions.R")

parametersDf <- expand.grid(rho = c(0.2,0.90),
                            p = c(2000),
                            SNR = c(0.2,1,2),
                            n = c(400), # this is the total train + test sample size
                            # nActive = c(300), # must be even because its being split among two modules
                            #n0 = 200,
                            cluster_distance = c("tom"),
                            Ecluster_distance = c("difftom"),
                            rhoOther = 0.6,
                            betaMean = c(1),
                            alphaMean = c(0.5,2),
                            betaE = 2,
                            includeInteraction = TRUE,
                            includeStability = TRUE,
                            distanceMethod = "euclidean",
                            clustMethod = "hclust",
                            #cutMethod = "gap",
                            cutMethod = "dynamic",
                            method = "average",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

parametersDf <- transform(parametersDf, n0 = n/2, nActive = p*0.05)
nSimScenarios <- nrow(parametersDf)


# MARS --------------------------------------------------------------------

# this loads an object called pp
load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/file19e17f13c434_2_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM.RData")
# this loads an object called pp2
load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/file19e150ebfb6e_2_PC_bouchard_sd_filter_10K_probes_corr_fisherScore.RData")

# to combine all the principal components
pcTrain_TOM <- pp$clustersAddon$PC
pp$clustersAll$nclusters
pp$clustersAddon$nclusters


dim(pcTrain_TOM)
head(pcTrain_TOM)
avgTrain_TOM <- pp$clustersAddon$averageExpr

# to combine all the principal components
pcTrain_corr <- pp2$clustersAddon$PC
dim(pcTrain_corr)
head(pcTrain_corr)
avgTrain_corr <- pp2$clustersAddon$averageExpr

varexp_PC1_TOM <- pp$clustersAddon$varExplained[seq(1, length(pp$clustersAddon$varExplained), by = 2)]
varexp_PC2_TOM <- pp$clustersAddon$varExplained[seq(2, length(pp$clustersAddon$varExplained), by = 2)]

varexp_PC1_corr <- pp2$clustersAddon$varExplained[seq(1, length(pp2$clustersAddon$varExplained), by = 2)]
varexp_PC2_corr <- pp2$clustersAddon$varExplained[seq(2, length(pp2$clustersAddon$varExplained), by = 2)]


dTOM <- data.frame(index = seq_len(length(varexp_PC1_TOM)), varexp_PC1_TOM, varexp_PC2_TOM) %>% 
  gather(type, value, -index) %>% 
  separate(type, c("measure", "PC", "type"))
dcorr <- data.frame(index = seq_len(length(varexp_PC1_corr)), varexp_PC1_corr, varexp_PC2_corr) %>% 
  gather(type, value, -index) %>% 
  separate(type, c("measure", "PC", "type"))

var_expl_data <- rbind(dTOM, dcorr)

p <- ggplot(var_expl_data, aes(x = index, y = value, color = PC))
p + geom_point(size=2) + facet_wrap(~type) + ylab("variance explained") + theme_bw()

max_heat <- max(c(max(pcTrain_TOM[which(pp$etrain==0),]),max(pcTrain_TOM[which(pp$etrain==1),])))
min_heat <- min(c(min(pcTrain_TOM[which(pp$etrain==0),]),min(pcTrain_TOM[which(pp$etrain==1),])))

library(ComplexHeatmap)
require(circlize)

cm <- colorRamp2(seq(min_heat, max_heat, length.out = 100), viridis(100))
ht1 = Heatmap(t(pcTrain_TOM[which(pp$etrain==0),]), 
              name = "E=0",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 0 : Age [4.8, 11.3]",
              # column_title = "Income_Level: 1-7",
              column_title = "NGD",
              show_row_names = FALSE)
ht2 = Heatmap(t(pcTrain_TOM[which(pp$etrain==1),]), 
              name = "E=1",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 1 : Age [11.3, 18]",
              # column_title = "Income_Level: 8-10",
              column_title = "GD",
              show_row_names = FALSE)
ht1 + ht2



kl <- prepare_data(data = cbind(pcTrain_TOM, pheno = Y[trainIndex], 
                                income = E[trainIndex]),
                   response = "pheno", exposure = "income")

kl$X %>% dim
kl$X[,1:155] %>% head


fit.earth <- earth::earth(x = kl$X[,1:154], y = kl$Y, pmethod = "forward", keepxy = TRUE, degree = 2, nfold = 10, trace = 4)
summary(fit.earth)
plot(fit.earth)
plotmo(fit.earth)
predict(fit.earth, newdata = kl$X)
evimp(fit.earth)


# resamples()
# 
# head(predict(marsTuned, kl$X))
# varImp(marsTuned)
# 
# earthPred <- predict(marsTuned, newdata = kl$X)
## The function ' postResample ' can be used to get the test set
## perforamnce values
# postResample(pred = earthPred, obs = kl$X)


fitControl <-  trainControl(method = "repeatedcv", 
                            number = 10, 
                            repeats = 3,
                            verboseIter = TRUE)

# Define the candidate models to test
marsGrid <- expand.grid(.degree = 1:2, .nprune = 2:20)
# Fix the seed so that the results can be reproduced
set.seed(1056)
marsTuned <- train(kl$X[,1:154], kl$Y,
                   method = "earth",
                   # Explicitly declare the candidate models to test
                   tuneGrid = marsGrid,
                   trControl = fitControl)
marsTuned$finalModel$bx
plot(marsTuned)
plotmo(marsTuned)

marsTuned$bestTune
fit.earth <- earth::earth(x = kl$X[,1:154], y = kl$Y, 
                          pmethod = "forward", keepxy = TRUE, degree = 2, nfold = 0, trace = 4)
plotmo(fit.earth)

set.seed(1056)
earthFit <- train(x = kl$X[,1:154], y = kl$Y,
                  method = "earth",
                  tuneLength = 10,
                  tuneGrid = marsGrid,
                  trControl = fitControl)

summary(earthFit)
plot(earthFit)
earthFit$results

set.seed(1056)
bagearthFit <- train(x = kl$X[,1:154], y = kl$Y,
                     method = "bagEarth",
                     tuneLength = 10,
                     trControl = fitControl)

plot(bagearthFit)

set.seed(1056)
lassoFit <- train(x = kl$X[,1:154], y = kl$Y,
                  method = "glmnet",
                  tuneLength = 10,
                  trControl = fitControl)

resamps <- resamples(list(MARS = earthFit,
                          bagMARS = bagearthFit,
                          GLMNET = lassoFit))

set.seed(1056)
relaxoFit <- train(x = kl$X[,1:154], y = kl$Y,
                   method = "relaxo",
                   tuneLength = 10,
                   trControl = fitControl)

set.seed(1056)
rbfSVM <- train(x = kl$X[,1:154], y = kl$Y,
                method = "svmRadial",
                tuneLength = 10,
                trControl = fitControl)

lassoFit
resamps <- resamples(list(MARS = earthFit,
                          bagMARS = bagearthFit,
                          GLMNET = lassoFit,
                          Relaxo = relaxoFit,
                          # Ridge = ridgeFit,
                          SVMrbf = rbfSVM))

bwplot(resamps)
summary(resamps)
dotplot(resamps)
splom(resamps)
modelDifferences <- diff(resamps)
?xyplot.resamples


predictors(bagearthFit)
predictors(earthFit)
predictors(lassoFit)
predictors(rbfSVM)

trellis.par.set(caretTheme())
plot(bagearthFit, type = c("g", "o"))


plot(earthFit)
summary(earthFit)
earthFit$finalModel
varImp(earthFit)
str(earthFit)

# or load pre-calculated results using:
load(url("http://caret.r-forge.r-project.org/exampleModels.RData"))

resamps <- resamples(list(CART = rpartFit,
                          CondInfTree = ctreeFit,
                          MARS = earthFit))

bwplot(resamps)
summary(resamps)
dotplot(resamps)






