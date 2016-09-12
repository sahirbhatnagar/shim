##################################
# R source code file for conducting NONLINEAR simulations
# on Mamouth cluster
# Git: this is on the eclust repo, sim2-modules-mammouth branch
# Created by Sahir,  Sept 10, 2016
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
                            p = c(1000),
                            SNR = c(0.2,1,2),
                            n = c(400), # this is the total train + test sample size
                            # nActive = c(300), # must be even because its being split among two modules
                            #n0 = 200,
                            cluster_distance = c("tom","corr"),
                            Ecluster_distance = c("difftom", "diffcorr"),
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
                            agglomerationMethod = "average",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

parametersDf <- transform(parametersDf, n0 = n/2, nActive = p*0.05)
parametersDf <- parametersDf[which(parametersDf$cluster_distance=="tom" & parametersDf$Ecluster_distance=="difftom" | 
                                     parametersDf$cluster_distance=="corr" & parametersDf$Ecluster_distance=="diffcorr"),]
nSimScenarios <- nrow(parametersDf)
parameterIndex <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))

parameterIndex = 5
simulationParameters <- parametersDf[parameterIndex,, drop = F]

print(simulationParameters)

## ---- generate-data ----
message("generating data")

p <- simulationParameters[,"p"];
n <- simulationParameters[,"n"];
n0 <- simulationParameters[,"n0"];
SNR <- simulationParameters[,"SNR"]
n1 <- n - n0
cluster_distance <- simulationParameters[,"cluster_distance"]
Ecluster_distance <- simulationParameters[,"Ecluster_distance"]
rhoOther <- simulationParameters[,"rhoOther"];
betaMean <- simulationParameters[,"betaMean"];
betaE <- simulationParameters[,"betaE"]
alphaMean <- simulationParameters[,"alphaMean"];
rho <- simulationParameters[,"rho"];
nActive <- simulationParameters[,"nActive"];
includeInteraction <- simulationParameters[,"includeInteraction"]
includeStability <- simulationParameters[,"includeStability"]
distanceMethod <- simulationParameters[,"distanceMethod"]
clustMethod <- simulationParameters[,"clustMethod"]
cutMethod <- simulationParameters[,"cutMethod"]
agglomerationMethod <- simulationParameters[,"agglomerationMethod"]
K.max <- simulationParameters[,"K.max"]
B <- simulationParameters[,"B"]

# in this simulation its blocks 3 and 4 that are important
# leaveOut:  optional specification of modules that should be left out 
# of the simulation, that is their genes will be simulated as unrelated 
# ("grey"). This can be useful when simulating several sets, in some which a module 
# is present while in others it is absent.
d0 <- simModule(n = n0, p = p, rho = c(0,0), exposed = FALSE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.01,
                maxCor = 1,
                corPower = 1,
                propNegativeCor = 0.3,
                backgroundNoise = 0.5,
                signed = FALSE,
                leaveOut = 1:4)

d1 <- simModule(n = n1, p = p, rho = c(rho, rho), exposed = TRUE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.4,
                maxCor = 1,
                corPower = 0.3,
                propNegativeCor = 0.3,
                backgroundNoise = 0.5,
                signed = FALSE)

# these should be the same. if they arent, its because I removed the red and
# green modules from the E=0 group
truemodule0 <- d0$setLabels
t0 <- table(truemodule0)
truemodule1 <- d1$setLabels
t1 <- table(truemodule1)
table(truemodule0,truemodule1)

# Convert labels to colors for plotting
moduleColors <- labels2colors(truemodule1)
table(moduleColors, truemodule1)

X <- rbind(d0$datExpr, d1$datExpr) %>%
  magrittr::set_colnames(paste0("Gene", 1:p)) %>%
  magrittr::set_rownames(paste0("Subject",1:n))

dim(X)
# pheatmap(cor(X))
# pheatmap(cor(d1$datExpr))
# pheatmap(cor(d0$datExpr))
# pheatmap(cor(d1$datExpr)-cor(d0$datExpr))
# pheatmap(WGCNA::TOMsimilarityFromExpr(X))
# pheatmap(WGCNA::TOMsimilarityFromExpr(d1$datExpr))
# pheatmap(WGCNA::TOMsimilarityFromExpr(d1$datExpr)-WGCNA::TOMsimilarityFromExpr(d0$datExpr))

# betaMainEffect <- vector("double", length = p)
# betaMainInteractions <- vector("double", length = p)
# # first assign random uniform to every gene in cluster 3 and 4,
# # then randomly remove so that thers only nActive left
# betaMainEffect[which(truemodule1 %in% 3:4)] <- runif(sum(truemodule1 %in% 3:4),
#                                                      betaMean - 0.1, betaMean + 0.1)
#
# # randomly set some coefficients to 0 so that there are only nActive non zero
# betaMainEffect <- replace(betaMainEffect,
#                           sample(which(truemodule1 %in% 3:4), sum(truemodule1 %in% 3:4) - nActive,
#                                  replace = FALSE), 0)
#
# betaMainInteractions[which(betaMainEffect!=0)] <- runif(nActive, alphaMean - 0.1, alphaMean + 0.1)
#
# beta <- c(betaMainEffect,
#           betaE,
#           betaMainInteractions)
# plot(beta)

betaMainEffect <- vector("double", length = p)
betaMainInteractions <- vector("double", length = p)

# the first nActive/2 in the 3rd block are active
betaMainEffect[which(truemodule1 %in% 3)[1:(nActive/2)]] <- runif(
  nActive/2, betaMean - 0.1, betaMean + 0.1)

# the first nActive/2 in the 4th block are active
betaMainEffect[which(truemodule1 %in% 4)[1:(nActive/2)]] <- runif(
  nActive/2, betaMean - 0.1, betaMean + 0.1)

betaMainInteractions[which(betaMainEffect!=0)] <- runif(nActive, alphaMean - 0.1, alphaMean + 0.1)

# must be in this order!!!! main effects, E, and then interactions... this order is being used
# by the generate_data function
beta <- c(betaMainEffect,
          betaE,
          betaMainInteractions)

#plot(beta)

result <- generate_data(p = p, X = X, 
                        beta = beta, include_interaction = includeInteraction,
                        cluster_distance = cluster_distance,
                        n = n, n0 = n0, 
                        eclust_distance = Ecluster_distance,
                        signal_to_noise_ratio = SNR,
                        distance_method = distanceMethod,
                        cluster_method = clustMethod,
                        cut_method = cutMethod,
                        agglomeration_method = agglomerationMethod,
                        K.max = K.max, B = B, nPC = 1)


# MARS --------------------------------------------------------------------

library(earth)
result %>% names


fit.earth <- earth::earth(x = result[["X_train"]], y = result[["Y_train"]], 
                          pmethod = "cv", keepxy = TRUE, 
                          degree = 2, trace = 4, nfold = 5)
# selected genes
coef(fit.earth)
dimnames(evimp(fit.earth))[[1]]
plot(fit.earth, which=1, col.rsq=0) # which=1 for Model Selection plot only (optional)



pp <- summary(fit.earth)

fit.earth %>% names
fit.earth$dirs



pp %>% names
pp$coefficients
plot(fit.earth)
plotmo(fit.earth)
predict(fit.earth, newdata = kl$X)




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






