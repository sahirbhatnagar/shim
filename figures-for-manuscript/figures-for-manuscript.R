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
                            agglomerationMethod = "average",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

parametersDf <- transform(parametersDf, n0 = n/2, nActive = p*0.05)
parametersDf <- parametersDf[which(parametersDf$cluster_distance=="tom" & parametersDf$Ecluster_distance=="difftom" | 
                                     parametersDf$cluster_distance=="corr" & parametersDf$Ecluster_distance=="diffcorr"),]
# parameterIndex <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))

parameterIndex = 2
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

#plot(beta)
X0 <- X[paste0("Subject",1:n0),]
X1 <- X[paste0("Subject",(n0+1):n),]

dim(X0)
dim(X1)

corrX <- cor(X)
class(corrX) <- c("similarity",class(corrX))

corrX0 <- cor(X0)
spearX0 <- WGCNA::cor(X0, method = "spearman")
class(corrX0) <- c("similarity",class(corrX0))

corrX1 <- cor(X1)
spearX1 <- WGCNA::cor(X1, method = "spearman")
class(corrX1) <- c("similarity",class(corrX1))

diffCorr <- abs(corrX1-corrX0)
class(diffCorr) <- c("similarity",class(diffCorr))

TOMX <- TOMsimilarityFromExpr(X)
class(TOMX) <- c("similarity",class(TOMX))
dimnames(TOMX)[[1]] <- dimnames(corrX)[[1]]
dimnames(TOMX)[[2]] <- dimnames(corrX)[[2]]

TOMX0 <- TOMsimilarityFromExpr(X0)
class(TOMX0) <- c("similarity",class(TOMX0))
dimnames(TOMX0)[[1]] <- dimnames(corrX)[[1]]
dimnames(TOMX0)[[2]] <- dimnames(corrX)[[2]]

TOMX1 <- TOMsimilarityFromExpr(X1)
class(TOMX1) <- c("similarity",class(TOMX1))
dimnames(TOMX1)[[1]] <- dimnames(corrX)[[1]]
dimnames(TOMX1)[[2]] <- dimnames(corrX)[[2]]


diffTOM <- abs(TOMX1 - TOMX0)
class(diffTOM) <- c("similarity",class(diffTOM))
dimnames(diffTOM)[[1]] <- dimnames(corrX)[[1]]
dimnames(diffTOM)[[2]] <- dimnames(corrX)[[2]]


# When alpha=2, the corr scor is proportional to  (p1+p2)/2 - p
alpha <- 2
Scorr <- abs(corrX0 + corrX1 - alpha * corrX)
class(Scorr) <- c("similarity", class(Scorr))

fisherScore <- fisherZ(n0 = n0,corrX0, n1 = n1, corrX1)

# on-off case as in Dettling 2005
Sonoff <- abs(spearX1 - spearX0)
class(Sonoff) <- c("similarity", class(Sonoff))


fig.w <- 1024
fig.h <- 768
fig.res <- 100

# path <- "/home/bhatnaga/coexpression/august2016simulation/linear/"
path <- ""

## ---- heat-corr-all ----
hc <- hclust(as.dist(1 - corrX), method = "average")
png(paste0(path, "figures-for-manuscript/corr_all.png"), width = fig.w, height = fig.h, res = fig.res)
plot(corrX, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc,
     active = as.numeric(betaMainEffect!=0))
dev.off()

# pheatmap(X, color = viridis(100))
# pheatmap(X0, color = viridis(100))
# pheatmap(X1, color = viridis(100))

## ---- heat-corr-e0 ----
hc <- hclust(as.dist(1 - corrX0), method = "average")
png(paste0(path, "figures-for-manuscript/corr_e0.png"), width = fig.w, height = fig.h, res = fig.res)
plot(corrX0, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc,
     active = as.numeric(betaMainEffect!=0))
dev.off()


## ---- heat-corr-e1 ----
hc <- hclust(as.dist(1 - corrX1), method = "average")
png(paste0(path, "figures-for-manuscript/corr_e1.png"), width = fig.w, height = fig.h, res = fig.res)
plot(corrX1, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc,
     active = as.numeric(betaMainEffect!=0))
dev.off()

# to reproduce the above you have use this code:
# plot(corrX1, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
#      clustering_distance_rows = as.dist(1 - corrX1),
#      clustering_distance_cols = as.dist(1 - corrX1))


## ---- heat-corr-diff ----
png(paste0(path, "figures-for-manuscript/corr_diff.png"), width = fig.w, height = fig.h, res = fig.res)
plot(diffCorr, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "average",
     clustering_distance_rows = dist(diffCorr),
     clustering_distance_cols = dist(diffCorr),
     active = as.numeric(betaMainEffect!=0))
dev.off()


## ---- heat-tom-all ----
hc <- hclust(as.dist(1 - TOMX), method = "average")
png(paste0(path, "figures-for-manuscript/tom_all.png"), width = fig.w, height = fig.h, res = fig.res)
plot(TOMX, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc,
     active = as.numeric(betaMainEffect!=0))
dev.off()

## ---- heat-tom-e0 ----
hc <- hclust(as.dist(1 - TOMX0), method = "average")
png(paste0(path, "figures-for-manuscript/tom_e0.png"), width = fig.w, height = fig.h, res = fig.res)
plot(TOMX0, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc,
     active = as.numeric(betaMainEffect!=0))
dev.off()

## ---- heat-tom-e1 ----
hc <- hclust(as.dist(1 - TOMX1), method = "average")
png(paste0(path, "figures-for-manuscript/tom_e1.png"), width = fig.w, height = fig.h, res = fig.res)
plot(TOMX1, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc,
     active = as.numeric(betaMainEffect!=0))
dev.off()


## ---- heat-tom-diff ----
hc <- hclust(as.dist(diffTOM), method = "average")
png(paste0(path, "figures-for-manuscript/tom_diff.png"), width = fig.w, height = fig.h, res = fig.res)
plot(diffTOM, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "average",
     clustering_distance_rows = dist(diffTOM),
     clustering_distance_cols = dist(diffTOM),
     active = as.numeric(betaMainEffect!=0))
dev.off()


## ---- fishers-zstat ----
png(paste0(path, "figures-for-manuscript/corr_fisher.png"), width = fig.w, height = fig.h, res = fig.res)
plot(fisherScore, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     active = as.numeric(betaMainEffect!=0))
dev.off()

## ---- cor-scor ----
png(paste0(path, "figures-for-manuscript/corr_scor.png"), width = fig.w, height = fig.h, res = fig.res)
plot(Scorr, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "average",
     active = as.numeric(betaMainEffect!=0))
dev.off()

## ---- on-off ----
png(paste0(path, "figures-for-manuscript/corr_onoff.png"), width = fig.w, height = fig.h, res = fig.res)
plot(Sonoff, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "average",
     active = as.numeric(betaMainEffect!=0))
dev.off()




