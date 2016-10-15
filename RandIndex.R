##################################
# R source code file for creating some of the WGCNA module plots 
# used for the manuscript
# based on the file found in:
# ~/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/bin/simulation/simDataWGCNA.R
# and ~/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/report/simulation-WGCNA.html
# Created by Sahir,  September 10th, 2016
# Updated:
# Notes:
# This is based on code from Network analysis book by Horvath
# it is being run on Mammouth
##################################

rm(list=ls())
source("packages.R")
source("functions.R")
options(digits = 4, scipen=999)

# source("/home/bhatnaga/coexpression/august2016simulation/linear/packages.R")
# source("/home/bhatnaga/coexpression/august2016simulation/linear/functions.R")

# source(paste(Sys.getenv("PBS_O_WORKDIR"),"packages.R", sep="/"))
# source(paste(Sys.getenv("PBS_O_WORKDIR"),"functions.R", sep="/"))

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
                            # agglomerationMethod = "average",
                            # agglomerationMethod = "ward.D2",
                            agglomerationMethod = "ward.D",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

parametersDf <- transform(parametersDf, n0 = n/2, nActive = p*0.05)
parametersDf <- parametersDf[which(parametersDf$cluster_distance=="tom" & parametersDf$Ecluster_distance=="difftom" | 
                                     parametersDf$cluster_distance=="corr" & parametersDf$Ecluster_distance=="diffcorr"),]
# parameterIndex <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))

for (parInd in c(3,12)) {

# parameterIndex = 3
  parameterIndex = parInd
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


message(sprintf("Rand Index for %s and true cluster membership is %g",
                cluster_distance,WGCNA::randIndex(table(result$clustersAll$cluster,truemodule1))))
message(sprintf("Rand Index for %s and true cluster membership is %g",
                Ecluster_distance,WGCNA::randIndex(table(result$clustersEclust$cluster,truemodule1))))
}
