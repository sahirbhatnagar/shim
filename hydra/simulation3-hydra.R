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

source(paste(Sys.getenv("PBS_O_WORKDIR"),"packages.R", sep="/"))
source(paste(Sys.getenv("PBS_O_WORKDIR"),"functions.R", sep="/"))

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
                            alphaMean = c(0.1),
                            betaE = 1,
                            includeInteraction = FALSE,
                            includeStability = TRUE,
                            distanceMethod = "euclidean",
                            clustMethod = "hclust",
                            #cutMethod = "gap",
                            cutMethod = "dynamic",
                            agglomerationMethod = "average",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

parametersDf <- transform(parametersDf, n0 = n/2, nActive = p*0.10)
parametersDf <- parametersDf[which(parametersDf$cluster_distance=="tom" & parametersDf$Ecluster_distance=="difftom" | 
                                     parametersDf$cluster_distance=="corr" & parametersDf$Ecluster_distance=="diffcorr"),]
nSimScenarios <- nrow(parametersDf)
parameterIndex <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))

# parameterIndex = 5
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

# this part is still necessary even though we arent actually using these 
# beta effects. This beta vector is passed to the generate_data_mars function
# in order to determine the active variables

betaMainEffect <- vector("double", length = p)

# the first nActive/2 in the 3rd block are active
betaMainEffect[which(truemodule1 %in% 3)[1:(nActive/2)]] <- runif(
  nActive/2, betaMean - 0.1, betaMean + 0.1)

# the first nActive/2 in the 4th block are active
betaMainEffect[which(truemodule1 %in% 4)[1:(nActive/2)]] <- runif(
  nActive/2, betaMean - 0.1, betaMean + 0.1)

# must be in this order!!!! main effects, E, and then interactions... this order is being used
# by the generate_data function
beta <- c(betaMainEffect,
          betaE)


result <- generate_data_mars(p = p, X = X, beta = beta, 
                             truemodule = truemodule1,
                             nActive = nActive,
                             include_interaction = includeInteraction,
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

# result %>% names
# 
# 
# 
# system.time(
# fit.earth <- earth::earth(x = result[["X_train"]], y = result[["Y_train"]], 
#                           keepxy = TRUE, nk = 1000,  
#                           degree = 2, trace = 4, nfold = 10, ncross = 3))
# # selected genes
# coef(fit.earth)
# dimnames(evimp(fit.earth))[[1]]
# get.used.pred.names(fit.earth)
# 
# plot(fit.earth, which=1, col.rsq=0) # which=1 for Model Selection plot only (optional)
# plot.earth.models(fit.earth$cv.list, which=1)
# plot(fit.earth)
# plot(fit.earth, which=1,
#      col.mean.infold.rsq="blue", col.infold.rsq="lightblue",
#      col.grsq=0, col.rsq=0, col.vline=0, col.oof.vline=0)
# 
# setdiff(get.used.pred.names(fit.earth),dimnames(evimp(fit.earth))[[1]])
# setdiff(dimnames(evimp(fit.earth))[[1]],get.used.pred.names(fit.earth))
# 
# dim(result[["X_test"]])
# dim(result[["DT_test"]])
# dtest <- result[["DT_test"]]
# colnames(dtest)[1] <- 'result[["Y_train"]]'
# 
# 
# yhat <- predict(fit.earth, newdata = result[["X_test"]])
# y <- result[["Y_test"]] # test response
# print(1 - sum((y - yhat)^2) / sum((y - mean(y))^2))
# sqrt(crossprod(y-yhat)/length(y))

print("starting MARS with interaction")
library(doMC)
registerDoMC(cores = 7)
mars_res <- mapply(mars_fun,
                  model = c("MARS"),
                  MoreArgs = list(x_train = result[["X_train"]],
                                  x_test = result[["X_test"]],
                                  y_train = result[["Y_train"]],
                                  y_test = result[["Y_test"]],
                                  stability = F,
                                  filter = F,
                                  filter_var = F,
                                  include_E = T,
                                  include_interaction = includeInteraction,
                                  s0 = result[["S0"]],
                                  p = p,
                                  true_beta = result[["beta_truth"]]),
                  SIMPLIFY = F,
                  USE.NAMES = F)

mars_res %>% unlist()

if (includeStability) {
  mars_stab <- mapply(function(model) mapply(mars_fun,
                                            x_train = result[["X_train_folds"]],
                                            y_train = result[["Y_train_folds"]],
                                            MoreArgs = list(stability = T,
                                                            model = model,
                                                            filter = F,
                                                            filter_var = F,
                                                            include_E = T,
                                                            include_interaction = includeInteraction,
                                                            s0 = result[["S0"]],
                                                            true_beta = result[["beta_truth"]],
                                                            p = p),
                                            SIMPLIFY = F),
                     #model = c("scad","mcp","lasso","elasticnet","ridge"),
                     model = c("MARS"),
                     SIMPLIFY = F,
                     USE.NAMES = F)
  
  
  # Make the combinations of list elements
  ll <- lapply(seq_along(mars_stab), function(i) combn(mars_stab[[i]], 2, simplify = F))
  
  
  mars_labels <- function(model) {
    paste0("mars_na_",model,"_yes_")
  }
  
  mars_labs <- mapply(mars_labels,
                     #model = c("scad","mcp","lasso","elasticnet","ridge"),
                     model = c("MARS"),
                     USE.NAMES = F)
  
  # Jaccard index
  mars_jacc <- lapply(seq_along(ll), function(j) {
    res <- mean(sapply(ll[[j]] , function(x) {
      A = x[[1]][coef.est != 0]$Gene
      B = x[[2]][coef.est != 0]$Gene
      if (length(A)==0 | length(B)==0) 0 else length(intersect(A,B))/length(union(A,B))
    }), na.rm = TRUE)
    names(res) <- paste0(mars_labs[[j]],"jacc")
    return(res)
  })
  
  mars_jacc %>% unlist
}

print("done MARS without interaction")




# Clust MARS --------------------------------------------------------------


# we will treat the clusters as fixed i.e., even if we filter, or
# do cross validation, the group labels are predetermined by the
# above clustering procedure
# This method is based on clusters derived without accounting for the
# environment

print("starting cluster and regress without interaction")

clust_res <- mapply(mars_clust_fun,
                    summary = rep(c("avg","pc"), each = 1),
                    model = c("MARS"),
                    MoreArgs = list(x_train = result[["X_train"]],
                                    x_test = result[["X_test"]],
                                    y_train = result[["Y_train"]],
                                    y_test = result[["Y_test"]],
                                    stability = F,
                                    filter = F,
                                    filter_var = F,
                                    include_E = T,
                                    include_interaction = includeInteraction,
                                    true_beta = result[["beta_truth"]],
                                    s0 = result[["S0"]],
                                    p = p,
                                    gene_groups = result[["clustersAll"]],
                                    clust_type = "clust",
                                    nPC = 1),
                    SIMPLIFY = F,
                    USE.NAMES = F)

# result %>% names
clust_res %>% unlist

if (includeStability) {
  clust_stab <- mapply(function(summary,
                                model) mapply(mars_clust_fun,
                                              x_train = result[["X_train_folds"]],
                                              y_train = result[["Y_train_folds"]],
                                              MoreArgs = list(stability = T,
                                                              x_test = result[["X_test"]],
                                                              summary = summary,
                                                              model = model,
                                                              filter = F,
                                                              filter_var = F,
                                                              include_E = T,
                                                              include_interaction = includeInteraction,
                                                              gene_groups = result[["clustersAll"]],
                                                              p = p,
                                                              true_beta = result[["beta_truth"]],
                                                              clust_type = "clust",
                                                              nPC = 1),
                                              SIMPLIFY = F),
                       summary = rep(c("avg","pc"), each = 1),
                       model = c("MARS"),
                       SIMPLIFY = F,
                       USE.NAMES = F)
  
  
  # Make the combinations of list elements
  ll <- lapply(seq_along(clust_stab), function(i) combn(clust_stab[[i]], 2, simplify = F))
  
  
  clust_labels <- function(summary, model) {
    paste0("clust",paste0("_",summary),paste0("_",model),"_","yes_")
  }
  
  clust_labs <- mapply(clust_labels,
                       summary = rep(c("avg","pc"), each = 1),
                       model = c("MARS"),
                       USE.NAMES = F)
  
  
  # Jaccard index
  clust_jacc <- lapply(seq_along(ll), function(j) {
    res <- mean(sapply(ll[[j]] , function(x) {
      A = x[[1]][coef.est != 0]$Gene
      B = x[[2]][coef.est != 0]$Gene
      if (length(A)==0 | length(B)==0) 0 else length(intersect(A,B))/length(union(A,B))
    }), na.rm = TRUE)
    names(res) <- paste0(clust_labs[[j]],"jacc")
    return(res)
  })
  
  clust_jacc %>% unlist
}

print("done clust and regress no interaction")


## ---- Ecluster-and-regress----

# we will treat the clusters as fixed i.e., even if we filter, or
# do cross validation, the group labels are predetermined by the
# above clustering procedure
# This method is based on clusters derived without accounting for the environment
# AND accoundting for the environment
# So we want to see if adding these clusters makes a difference from what people
# might usually do, i.e just based on correlation without E

message("starting Environment cluster and regress without interaction")

Eclust_res <- mapply(mars_clust_fun,
                     summary = rep(c("avg","pc"), each = 1),
                     model = c("MARS"),
                     MoreArgs = list(x_train = result[["X_train"]],
                                     x_test = result[["X_test"]],
                                     y_train = result[["Y_train"]],
                                     y_test = result[["Y_test"]],
                                     stability = F,
                                     filter = F,
                                     filter_var = F,
                                     include_E = T,
                                     include_interaction = includeInteraction,
                                     s0 = result[["S0"]],
                                     p = p,
                                     true_beta = result[["beta_truth"]],
                                     gene_groups = result[["clustersAddon"]],
                                     clust_type = "Eclust",
                                     nPC = 1),
                     SIMPLIFY = F,
                     USE.NAMES = F)

# Eclust_res %>% names
# options(digits = 2, scipen=999)
mars_res %>% unlist()
clust_res %>% unlist()
Eclust_res %>% unlist

if (includeStability) {
  Eclust_stab <- mapply(function(summary,
                                 model) mapply(mars_clust_fun,
                                               x_train = result[["X_train_folds"]],
                                               y_train = result[["Y_train_folds"]],
                                               MoreArgs = list(stability = T,
                                                               x_test = result[["X_test"]],
                                                               summary = summary,
                                                               model = model,
                                                               filter = F,
                                                               filter_var = F,
                                                               include_E = T,
                                                               include_interaction = includeInteraction,
                                                               gene_groups = result[["clustersAddon"]],
                                                               p = p,
                                                               true_beta = result[["beta_truth"]],
                                                               clust_type = "Eclust",
                                                               nPC = 1),
                                               SIMPLIFY = F),
                        summary = rep(c("avg","pc"), each = 1),
                        model = c("MARS"),
                        SIMPLIFY = F,
                        USE.NAMES = F)
  
  
  # Make the combinations of list elements
  ll <- lapply(seq_along(Eclust_stab), function(i) combn(Eclust_stab[[i]], 2, simplify = F))
  
  Eclust_labels <- function(summary, model) {
    paste0("Eclust",paste0("_",summary),paste0("_",model),"_","yes_")
  }
  
  Eclust_labs <- mapply(Eclust_labels,
                        summary = rep(c("avg","pc"), each = 1),
                        model = c("MARS"),
                        USE.NAMES = F)
  
  
  # Jaccard index
  Eclust_jacc <- lapply(seq_along(ll), function(j) {
    res <- mean(sapply(ll[[j]] , function(x) {
      A = x[[1]][coef.est != 0]$Gene
      B = x[[2]][coef.est != 0]$Gene
      if (length(A)==0 | length(B)==0) 0 else length(intersect(A,B))/length(union(A,B))
    }), na.rm = TRUE)
    names(res) <- paste0(Eclust_labs[[j]],"jacc")
    return(res)
  })
  
  Eclust_jacc %>% unlist
}

print("done Environment clust and regress with interaction")



## ---- final-results ----

final_results <- if (includeStability) {
  c(simulationParameters,
    # uni_res,uni_mean_stab, "uni_jacc" = uni_jacc,
    clust_res, clust_jacc,
    mars_res, mars_jacc, 
    Eclust_res, Eclust_jacc) %>% unlist } else {
      c(simulationParameters, 
        # uni_res,
        clust_res,
        mars_res,
        Eclust_res)  %>% unlist
    }

final_results %>% t %>% as.data.frame()


filename <- tempfile(pattern = paste0(sprintf("%s_%s_rho%.2f_p%1.0f_SNR%.2f_n%1.0f_s0%1.0f_beta%.2f_alpha%.2f",
                                              cluster_distance,Ecluster_distance,rho,p,SNR, n, nActive, betaMean, alphaMean),"_"),
                     tmpdir = paste(Sys.getenv("PBS_O_WORKDIR"),"results", sep="/"))
#tmpdir = "/home/bhatnaga/coexpression/may2016simulation/sim2-modules-mammouth/results/")
write.table(final_results %>% t %>% as.data.frame(),
            file = filename,
            quote = F,
            row.names = F,
            col.names = F)

write.table(final_results %>% t %>% as.data.frame() %>% colnames(),
            file = paste(Sys.getenv("PBS_O_WORKDIR"),"colnames_stab_hydra-sim3-sept27.txt",sep="/"),
            # file  = filename,
            quote = F,
            row.names = F, col.names = F)

# final_results %>% t %>% as.data.frame() %>% colnames()


# # Mars simulation ESL page 327 ---------------------------------------------------------
# 
# result[["X_train"]] %>% dim
# 
# 
# # Scenario 1
# 
# X <- replicate(3, rnorm(100))
# dimnames(X)[[2]] <- c("X1","X2","e")
# Y <- pmax(X[,"X1"]-1,0) + pmax(X[,"X1"]-1,0) * pmax(X[,"X2"]-0.8,0) + 0.12 * X[,"e"]
# DT <- cbind(Y,X)
# 
# fit1 <- earth::earth(x = X[,c("X1","X2")], 
#                      y = Y,
#                      keepxy = TRUE,
#                      pmethod = "forward",
#                      # nk = 1000,
#                      degree = 2, 
#                      trace = 4,
#                      nfold = 10, ncross = 5)
# summary(fit1)
# plotmo(fit1, which = 3)
# plot(fit1, which=1, col.rsq=0) # which=1 for Model Selection plot only (optional)
# plot.earth.models(fit1$cv.list, which=1)
# plot(fit1)
# plot(fit1, which=1,
#      col.mean.infold.rsq="blue", col.infold.rsq="lightblue",
#      col.grsq=0, col.rsq=0, col.vline=0, col.oof.vline=0)
# 
# 
# # Figure Used in manuscript ------------------------------------------------------
# 
# hingeprod <- function(x1, x2) {
#   0.1*(x1 + x2+1) + 4 * pmax(x1-0.01, 0) * pmax(x2-0.05, 0)
#   # 10*sin(pi * x*y)
#   # tan(pmax(x,0) * pmax(y,0))
#   # 0.1*exp(4*x) + 4/(1+exp(-20*(y-0.5)))
# }
# 
# 
# linprod <- function(x1, x2) {
#   0.1*(x1 + x2)
#   # 10*sin(pi * x*y)
#   # tan(pmax(x,0) * pmax(y,0))
#   # 0.1*exp(4*x) + 4/(1+exp(-20*(y-0.5)))
# }
# 
# outseq <- function(x, y, fun){
#   x1r <- range(x)
#   x1seq <- seq(x1r[1], x1r[2], length.out = 25)
#   x2r <- range(y)
#   x2seq <- seq(x2r[1], x2r[2], length.out = 25)
#   return(list(x = x1seq, y = x2seq, z = outer(x1seq, x2seq, fun)))
# }
#
# result[["S0"]]
# x1 <- result[["X_train"]][, result[["S0"]][1:25]]
# dim(x1)
# u1 <- svd(x1)$u[,1]
# # u1 <- apply(x1, 1, mean)
#
# x2 <- result[["X_train"]][, result[["S0"]][26:50]]
# dim(x2)
# u2 <- svd(x2)$u[,1]
# # u2 <- apply(x2, 1, mean)
#
# reslin <- outseq(u1,u2, linprod)
# reshinge <- outseq(u1,u2, hingeprod)
# persp(res[["y"]], res[["x"]], res[["z"]])
# persp(res[["y"]], res[["x"]], res[["z"]],
#       theta = 10, phi = 30, expand = 0.8, col = "lightblue", ltheta = 120, shade = .01, ticktype = "simple",    xlab = "X1", ylab = "X2", zlab = "(X1-1)_+ * (X2-0.8)_+")

# savepdf("~/git_repositories/eclust-simulation-aug2016/sim3-persp.pdf")
# par(mfrow=c(1,2))
# drape.plot(reslin[["y"]], reslin[["x"]], reslin[["z"]], col=colorRampPalette(brewer.pal(9, "Reds"))(100),
#            horizontal = FALSE, add.legend = F, zlim = range(reshinge[["z"]]), zlim2 = range(reshinge[["z"]]),
#            theta = 50, phi = 30, expand = 0.75, ltheta = 120,
#            xlab = "U1", ylab = "U2", zlab = "Y")#, main = "E = 0")
# drape.plot(reshinge[["y"]], reshinge[["x"]], reshinge[["z"]], col=colorRampPalette(brewer.pal(9, "Reds"))(100),
#            horizontal = FALSE, add.legend = F, zlim = range(reshinge[["z"]]), zlim2 = range(reshinge[["z"]]),
#            theta = 50, expand = 0.75, phi = 30, ltheta = 120,
#            xlab = "U1", ylab = "U2", zlab = "Y")#, main = "E = 1")
# dev.off()



hingeprod <- function(x1, x2) {
  1*(x1 + 1) + 1 * x2
  # 10*sin(pi * x*y)
  # tan(pmax(x,0) * pmax(y,0))
  # 0.1*exp(4*x) + 4/(1+exp(-20*(y-0.5)))
}


linprod <- function(x1, x2) {
  1*(x1 + 1)
  # 10*sin(pi * x*y)
  # tan(pmax(x,0) * pmax(y,0))
  # 0.1*exp(4*x) + 4/(1+exp(-20*(y-0.5)))
}

outseq <- function(x, y, fun){
  x1r <- range(x)
  x1seq <- seq(x1r[1], x1r[2], length.out = 25)
  x2r <- range(y)
  x2seq <- seq(x2r[1], x2r[2], length.out = 25)
  return(list(x = x1seq, y = x2seq, z = outer(x1seq, x2seq, fun)))
}


reslin <- outseq(x_main,x_inter, linprod)
reshinge <- outseq(x_main,x_inter, hingeprod)

savepdf("~/git_repositories/eclust-simulation-aug2016/sim3-persp-Qi.pdf")
par(mfrow=c(1,2))
drape.plot(reslin[["y"]], reslin[["x"]], reslin[["z"]], col=colorRampPalette(brewer.pal(9, "Reds"))(100),
           horizontal = FALSE, add.legend = F, zlim = range(reshinge[["z"]]), zlim2 = range(reshinge[["z"]]),
           theta = 50, phi = 30, expand = 0.75, ltheta = 120,
           xlab = "1st PC", ylab = "f(Q)", zlab = "Y")#, main = "E = 0")
drape.plot(reshinge[["y"]], reshinge[["x"]], reshinge[["z"]], col=colorRampPalette(brewer.pal(9, "Reds"))(100),
           horizontal = FALSE, add.legend = F, zlim = range(reshinge[["z"]]), zlim2 = range(reshinge[["z"]]),
           theta = 50, expand = 0.75, phi = 30, ltheta = 120,
           xlab = "1st PC", ylab = "f(Q)", zlab = "Y")#, main = "E = 1")
dev.off()


