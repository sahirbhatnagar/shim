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

source("packages.R")
source("functions.R")

parametersDf <- expand.grid(rho = c(0.1,0.35,0.75,0.95),
                            blocksize = 50,
                            p = c(500, 1000, 3000),
                            n = 400, n0 = 200,
                            cluster = "diffcorr",
                            rhoOther = 0.6,
                            betaMean = 4,
                            betaE = 5,
                            alphaMean = 2,
                            nActive = 20,
                            includeInteraction = TRUE,
                            includeStability = TRUE, stringsAsFactors = FALSE)

#str(parametersDf)
parametersDf <- transform(parametersDf, nBlocks = p/blocksize)
parameterIndex <- commandArgs(trailingOnly = T)

parameterIndex = 7
simulationParameters <- parametersDf[parameterIndex,, drop = F]

## ---- generate-data ----
message("generating data")

p <- simulationParameters[,"p"];
n <- simulationParameters[,"n"];
n0 <- simulationParameters[,"n0"];
n1 <- n - n0
eclustMatrix <- simulationParameters[,"cluster"]
rhoOther <- simulationParameters[,"rhoOther"];
blocksize <- simulationParameters[,"blocksize"];
betaMean <- simulationParameters[,"betaMean"];
betaE <- simulationParameters[,"betaE"]
alphaMean <- simulationParameters[,"alphaMean"];
nBlocks <- simulationParameters[,"nBlocks"];
rho <- simulationParameters[,"rho"];
nActive <- simulationParameters[,"nActive"];
includeInteraction <- simulationParameters[,"includeInteraction"]
includeStability <- simulationParameters[,"includeStability"]

n0 = 30
n1 = 70
n = 100
p = 500
rho = 0.75

d0 <- simModule(n0, p, c(rho,-rho), exposed = FALSE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.3,
                maxCor = 1,
                corPower = 0.3,
                propNegativeCor = 0.1,
                backgroundNoise = 0.2,
                signed = TRUE,
                leaveOut = 3:4)

d1 <- simModule(n1, p, c(rho, rho), exposed = TRUE,
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
table(truemodule0,truemodule1)

# Convert labels to colors for plotting
moduleColors <- labels2colors(truemodule1)

X <- rbind(d0$datExpr, d1$datExpr) %>%
  magrittr::set_colnames(paste0("Gene", 1:p)) %>%
  magrittr::set_rownames(paste0("Subject",1:n))
