##################################
# R source code file for conducting simulations
# on HYDRA cluster
# Git: this is on the eclust repo, simulation branch
# Created by Sahir,  April 2, 2016
# Updated:
# Notes:
# This is a modified simulation and different from simulation1
# Its based on code from Network analysis book by Horvath
##################################

rm(list=ls())
source("packages.R")
source("functions.R")

parametersDf <- expand.grid(rho = c(0.1,0.35,0.75,0.95),
                            p = c(100, 500, 1000, 3000),
                            n = 500, n0 = 250,
                            cluster_distance = c("corr"),
                            rhoOther = 0.6,
                            betaMean = 4,
                            betaE = 5,
                            alphaMean = 2,
                            nActive = 50,
                            includeInteraction = TRUE,
                            includeStability = TRUE,
                            distanceMethod = "euclidean",
                            clustMethod = "hclust",
                            cutMethod = "gap",
                            method = "complete",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

#str(parametersDf)
parametersDf <- transform(parametersDf, nBlocks = p/blocksize)
parameterIndex <- commandArgs(trailingOnly = T)

parameterIndex = 11
simulationParameters <- parametersDf[parameterIndex,, drop = F]

## ---- generate-data ----
message("generating data")

p <- simulationParameters[,"p"];
n <- simulationParameters[,"n"];
n0 <- simulationParameters[,"n0"];
n1 <- n - n0
cluster_distance <- simulationParameters[,"cluster_distance"]
rhoOther <- simulationParameters[,"rhoOther"];
betaMean <- simulationParameters[,"betaMean"];
betaE <- simulationParameters[,"betaE"]
alphaMean <- simulationParameters[,"alphaMean"];
nBlocks <- simulationParameters[,"nBlocks"];
rho <- simulationParameters[,"rho"];
nActive <- simulationParameters[,"nActive"];
includeInteraction <- simulationParameters[,"includeInteraction"]
includeStability <- simulationParameters[,"includeStability"]
distanceMethod <- simulationParameters[,"distanceMethod"]
clustMethod <- simulationParameters[,"clustMethod"]
cutMethod <- simulationParameters[,"cutMethod"]
method <- simulationParameters[,"method"]
K.max <- simulationParameters[,"K.max"]
B <- simulationParameters[,"B"]

# in this simulation its blocks 3 and 4 that are important
d0 <- simModule(n = n0, p = p, rho = c(rho,-rho), exposed = FALSE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.3,
                maxCor = 1,
                corPower = 0.3,
                propNegativeCor = 0.1,
                backgroundNoise = 0.2,
                signed = TRUE,
                leaveOut = 3:4)

d1 <- simModule(n = n1, p = p, rho = c(rho, rho), exposed = TRUE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.3,
                maxCor = 1,
                corPower = 0.3,
                propNegativeCor = 0.1,
                backgroundNoise = 0.2,
                signed = TRUE)

# these should be the same. if they arent, its because I removed the red and
# green modules from the E=0 group
truemodule0 <- d0$setLabels
t0 <- table(truemodule0)
truemodule1 <- d1$setLabels
t1 <- table(truemodule1)
table(truemodule0,truemodule1)

# Convert labels to colors for plotting
moduleColors <- labels2colors(truemodule1)

X <- rbind(d0$datExpr, d1$datExpr) %>%
  magrittr::set_colnames(paste0("Gene", 1:p)) %>%
  magrittr::set_rownames(paste0("Subject",1:n))

dim(X)


betaMainEffect <- vector("double", length = p)
betaMainInteractions <- vector("double", length = p)
# first assign random uniform to every gene in cluster 3 and 4,
# then randomly remove so that thers only nActive left
betaMainEffect[which(truemodule1 %in% 3:4)] <- runif(sum(truemodule1 %in% 3:4),
                                           betaMean - 0.1, betaMean + 0.1)

# randomly set some coefficients to 0 so that there are only nActive non zero
betaMainEffect <- replace(betaMainEffect,
                          sample(which(truemodule1 %in% 3:4), sum(truemodule1 %in% 3:4) - nActive,
                                 replace = FALSE), 0)

betaMainInteractions[which(betaMainEffect!=0)] <- runif(nActive, alphaMean - 0.1, alphaMean + 0.1)

beta <- c(betaMainEffect,
          betaE,
          betaMainInteractions)
plot(beta)

result <- generate_data(p = p, n = n, n0 = n0, X = X,
                        beta = beta, include_interaction = includeInteraction,
                        cluster_distance = cluster_distance,
                        EclustAddDistance = "fisherScore",
                        signal_to_noise_ratio = 1/2,
                        distanceMethod = distanceMethod,
                        clustMethod = clustMethod,
                        cutMethod = cutMethod,
                        method = method,
                        K.max = K.max, B = B)









