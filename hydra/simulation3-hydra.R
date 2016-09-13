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



generate_data_mars <- function(p, X, beta,
                          cluster_distance = c("corr", "corr0", "corr1", "tom",
                                               "tom0", "tom1", "diffcorr",
                                               "difftom","corScor", "tomScor",
                                               "fisherScore"),
                          n, n0, include_interaction = F,
                          signal_to_noise_ratio = 1,
                          eclust_distance = c("fisherScore", "corScor", "diffcorr",
                                              "difftom"),
                          cluster_method = c("hclust", "protoclust"),
                          cut_method = c("dynamic","gap", "fixed"),
                          distance_method = c("euclidean","maximum", "manhattan",
                                              "canberra", "binary", "minkowski"),
                          n_clusters,
                          agglomeration_method = c("complete", "average", "ward.D2",
                                                   "single", "ward.D", "mcquitty",
                                                   "median", "centroid"),
                          nPC = 1, 
                          K.max = 10, B = 10) {
  
  # p = p; X = X ; beta = beta
  # n = n; n0 = n0
  # cluster_distance = "corr"
  # include_interaction = F
  # signal_to_noise_ratio = 0.5
  # cluster_method = "hclust" ; cut_method = "dynamic";agglomeration_method="complete";
  # distance_method = "euclidean"
  # eclust_distance = "diffcorr"; nPC = 1
  
  
  agglomeration_method <- match.arg(agglomeration_method)
  cut_method <- match.arg(cut_method)
  cluster_method <- match.arg(cluster_method)
  distance_method <- match.arg(distance_method)
  cluster_distance <- match.arg(cluster_distance)
  eclust_distance <- match.arg(eclust_distance)
  
  
  names(beta) <- if (include_interaction) {
    c(paste0("Gene",1:p),"E", paste0("Gene",1:p,":E"))
  } else c(paste0("Gene",1:p),"E")
  
  # total true beta vector: this includes all the betas for the genes, then the
  # environment beta, then their interactions if interaction is true.
  # This is used to calculate the model error. This is the same as beta,
  # but in matrix form
  beta_truth <- as.matrix(beta)
  
  # Gene names belonging to the active set
  S0 <- names(beta)[which(beta != 0)]
  
  n1 <- n - n0
  
  message("Creating data and simulating response")
  
  DT <- as.data.frame(sim_data(n = n, n0 = n0, p = p, genes = X,
                               include_interaction = include_interaction,
                               E = c(rep(0,n0), rep(1, n1)),
                               beta = beta,
                               signal_to_noise_ratio = signal_to_noise_ratio))
  dim(DT)
  
  Y <- as.matrix(DT[,"Y"])
  
  #remove response from X0 and X1
  X0 <- as.matrix(DT[which(DT$E == 0),-1])
  X1 <- as.matrix(DT[which(DT$E == 1),-1])

  # partition-data
  trainIndex <- caret::createDataPartition(DT$E, p = .5, list = FALSE, times = 1)
  DT_train <- DT[trainIndex,]
  DT_test <- DT[-trainIndex,]
  
  # X_train and X_test contain the environment variable
  X_train <- DT_train[,-1] %>% as.matrix
  Y_train <- DT_train[, 1]
  X_test <- DT_test[,-1] %>% as.matrix
  Y_test <- DT_test[, 1]
  
  mse_null <- crossprod(mean(Y_test) - Y_test)/length(Y_test)
  
  # gene expression data
  genes_e0 <- DT_train[which(DT_train$E == 0),paste0("Gene",1:p)] %>% as.matrix
  genes_e1 <- DT_train[which(DT_train$E == 1),paste0("Gene",1:p)] %>% as.matrix
  genes_all <- rbind(genes_e0,genes_e1)
  
  message("Calculating similarity matrices")
  
  # gene expression data
  genes_all_test <- DT_test[,paste0("Gene",1:p)] %>% as.matrix
  
  corr_train_e0 <- WGCNA::cor(genes_e0)
  corr_train_e1 <- WGCNA::cor(genes_e1)
  corr_train_diff <- abs(corr_train_e1 - corr_train_e0)
  corr_train_all <- WGCNA::cor(genes_all)
  
  tom_train_e0 <- WGCNA::TOMsimilarityFromExpr(genes_e0)
  dimnames(tom_train_e0)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_e0)[[2]] <- dimnames(corr_train_all)[[2]]
  
  tom_train_e1 <- WGCNA::TOMsimilarityFromExpr(genes_e1)
  dimnames(tom_train_e1)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_e1)[[2]] <- dimnames(corr_train_all)[[2]]
  
  tom_train_diff <- abs(tom_train_e1 - tom_train_e0)
  dimnames(tom_train_diff)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_diff)[[2]] <- dimnames(corr_train_all)[[2]]
  
  tom_train_all <- WGCNA::TOMsimilarityFromExpr(genes_all)
  dimnames(tom_train_all)[[1]] <- dimnames(corr_train_all)[[1]]
  dimnames(tom_train_all)[[2]] <- dimnames(corr_train_all)[[2]]
  
  
  
  
  # corScor and Fisher Score matrices
  alpha <- 2
  Scorr <- abs(corr_train_e0 + corr_train_e1 - alpha * corr_train_all)
  class(Scorr) <- c("similarity", class(Scorr))
  
  # Stom <- abs(tom_train_e1 + tom_train_e0 - alpha * tom_train_all)
  # class(Stom) <- c("similarity", class(Stom))
  
  fisherScore <- fisherZ(n0 = n0, cor0 = corr_train_e0,
                         n1 = n1, cor1 = corr_train_e1)
  
  # class(tom_train_all) <- append(class(tom_train_all), "similarity")
  # class(tom_train_diff) <- append(class(tom_train_diff), "similarity")
  # class(tom_train_e1) <- append(class(tom_train_e1), "similarity")
  # class(tom_train_e0) <- append(class(tom_train_e0), "similarity")
  class(corr_train_all) <- append(class(corr_train_all), "similarity")
  class(corr_train_diff) <- append(class(corr_train_diff), "similarity")
  class(corr_train_e1) <- append(class(corr_train_e1), "similarity")
  class(corr_train_e0) <- append(class(corr_train_e0), "similarity")
  
  message("Creating CV folds from training data")
  
  # Folds for Cross validation
  folds_train <- caret::createFolds(Y_train, k = 10, list = T)
  DT_train_folds <- lapply(folds_train, function(i) DT_train[-i,])
  X_train_folds <- lapply(DT_train_folds, function(i) i[,-grep("Y",colnames(i))])
  Y_train_folds <- lapply(DT_train_folds, function(i) i[,grep("Y",colnames(i))])
  
  message(sprintf("Calculating number of clusters based on %s using %s with %s
                  linkage and the %s to determine the number of clusters",
                  cluster_distance, cluster_method, agglomeration_method, cut_method))
  
  # clusters based on cluster_distance argument
  similarity <- switch(cluster_distance,
                       corr = corr_train_all,
                       corr0 = corr_train_e0,
                       corr1 = corr_train_e1,
                       diffcorr = corr_train_diff,
                       difftom = tom_train_diff,
                       tom0 = tom_train_e0,
                       tom1 = tom_train_e1,
                       tom = tom_train_all,
                       corScor = Scorr,
                       tomScor = Stom,
                       fisherScore = fisherScore)
  
  # results for clustering, PCs and averages for each block
  # the only difference here is the distance_method arg
  res <- if (cluster_distance %in% c("diffcorr","difftom",
                                     "corScor", "tomScor","fisherScore")) {
    clusterSimilarity(x = similarity,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      distanceMethod = distance_method,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  } else {
    clusterSimilarity(x = similarity,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  }
  
  message(paste("Calculating number of environment clusters based on ",
                eclust_distance))
  
  # clusters based on eclust_distance
  similarityEclust <- switch(eclust_distance,
                             corr = corr_train_all,
                             corr0 = corr_train_e0,
                             corr1 = corr_train_e1,
                             diffcorr = corr_train_diff,
                             difftom = tom_train_diff,
                             tom0 = tom_train_e0,
                             tom1 = tom_train_e1,
                             tom = tom_train_all,
                             corScor = Scorr,
                             tomScor = Stom,
                             fisherScore = fisherScore)
  
  
  resEclust <- if (eclust_distance %in% c("diffcorr","difftom",
                                            "corScor", "tomScor","fisherScore")) {
    clusterSimilarity(x = similarityEclust,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      distanceMethod = distance_method,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  } else {
    clusterSimilarity(x = similarityEclust,
                      expr = genes_all,
                      exprTest = genes_all_test,
                      clustMethod = cluster_method,
                      cutMethod = cut_method,
                      method = agglomeration_method,
                      K.max = K.max, B = B, nClusters = nClusters, nPC = nPC)
  }
  
  
  # we need to combine the cluster information here
  # this is based on cluster_distance only
  clustersAll <- copy(res$clusters)
  n_clusters_All <- res$pcInfo$nclusters
  
  message(sprintf("There are %d clusters derived from the %s similarity matrix",
                  n_clusters_All, cluster_distance))
  
  # this is based on eclust_distance only
  n_clusters_Eclust <- resEclust$pcInfo$nclusters
  clustersEclust <- copy(resEclust$clusters)
  
  message(sprintf("There are %d clusters derived from the %s environment similarity matrix",
                  n_clusters_Eclust, eclust_distance))
  
  # this is based on both
  n_clusters_Addon <- n_clusters_All + n_clusters_Eclust
  
  message(sprintf("There are a total of %d clusters derived from the %s
                  similarity matrix and the %s environment similarity matrix",
                  n_clusters_Addon,cluster_distance,eclust_distance))
  
  # check if any of the cluster numbers in clustersEclust are 0
  # if there are, then add n_clusters+1 to each module number in
  # clustersEclust, else just add n_clusters. this is to control for the
  # situation where there are some clusters numbers of 0 which would cause
  # identical cluster numbers in the clusters and clustersEclust data
  if (clustersEclust[,any(cluster==0)]) {
    clustersEclust[,cluster := cluster + n_clusters_All + 1 ]
  } else {
    clustersEclust[,cluster := cluster + n_clusters_All ]
  }
  
  # this contains the clusters from the cluster_distance (e.g. corr matrix)
  # and the clusters from the eclust_distance (e.g. fisherScore)
  clustersAddon <- rbindlist(list(clustersAll, clustersEclust))
  
  # need to calculate penalty factors for group lasso
  # I put all main effects and interactions of a given module in the same group
  # and the size of the penalty factor is sqrt(size of module), where the
  # size of the module includes both main and interaction effects
  # environment should get penalized, in the original simulation 1
  # it was not being penalized which is maybe why it was performing well
  if (include_interaction) {
    
    gene_groups = copy(clustersAll)
    gene_groups[, gene := paste0(gene,":E")]
    gene_groups <- rbind(clustersAll,gene_groups) %>% setkey(cluster)
    
    pf_temp <- gene_groups[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)
    
    gene_groups_inter <- rbind(pf_temp[gene_groups],
                               data.table(cluster = n_clusters_All, N = 1,
                                          pf = 1, gene = "E", module = "empty"))
    # gglasso needs groups number consecutively 1, 2,3 ...
    gene_groups_inter[, cluster:=cluster+1]
    setkey(gene_groups_inter, cluster)
    
    gene_groups_Addon = copy(clustersAddon)
    gene_groups_Addon[, gene := paste0(gene,":E")]
    gene_groups_Addon <- rbind(clustersAddon, gene_groups_Addon) %>% setkey(cluster)
    
    pf_temp_Addon <- gene_groups_Addon[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)
    
    gene_groups_inter_Addon <- rbind(pf_temp_Addon[gene_groups_Addon],
                                     data.table(cluster = n_clusters_Addon, N = 1,
                                                pf = 1, gene = "E", module = "empty"))
    # gglasso needs groups number consecutively 1, 2,3 ...
    gene_groups_inter_Addon[, cluster:=cluster+1]
    setkey(gene_groups_inter_Addon, cluster)
  }
  
  DT <- DT %>% as.matrix
  class(DT) <- append(class(DT),"eset")
  
  result <- list(beta_truth = beta_truth,
                 similarity = similarity,
                 similarityEclust = similarityEclust,
                 DT = DT,
                 Y = Y, X0 = X0, X1 = X1, X_train = X_train, X_test = X_test,
                 Y_train = Y_train, Y_test = Y_test, DT_train = DT_train,
                 DT_test = DT_test, S0 = S0,
                 n_clusters_All = n_clusters_All,
                 n_clusters_Eclust = n_clusters_Eclust,
                 n_clusters_Addon = n_clusters_Addon,
                 clustersAll = clustersAll,
                 clustersAddon = clustersAddon,
                 clustersEclust = clustersEclust,
                 gene_groups_inter = if (include_interaction) gene_groups_inter else NULL,
                 gene_groups_inter_Addon = if (include_interaction) gene_groups_inter_Addon else NULL,
                 tom_train_all = tom_train_all, tom_train_diff = tom_train_diff,
                 tom_train_e1 = tom_train_e1,tom_train_e0 = tom_train_e0,
                 corr_train_all = corr_train_all,
                 corr_train_diff = corr_train_diff,
                 corr_train_e1 = corr_train_e1, corr_train_e0 = corr_train_e0,
                 fisherScore = fisherScore,
                 corScor = Scorr,
                 # corTom = Stom,
                 mse_null = mse_null, DT_train_folds = DT_train_folds,
                 X_train_folds = X_train_folds, Y_train_folds = Y_train_folds)
  return(result)
}



mars_fun <- function(x_train,
                     x_test,
                     y_train,
                     y_test,
                     s0,
                     model = c("MARS"),
                     true_beta,
                     topgenes = NULL,
                     stability = F,
                     filter = F,
                     include_E = F,
                     include_interaction = F,
                     p = 1000,
                     filter_var = F, ...){
  
  # stability = F; x_train = result[["X_train"]] ; x_test = result[["X_test"]] ;
  # y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = T;
  # s0 = result[["S0"]]; p = p ;
  # model = "MARS"; topgenes = NULL; true_beta = result[["beta_truth"]]
  
  # stability = F; x_train = result_interaction[["X_train"]] ; x_test = result_interaction[["X_test"]] ;
  # y_train = result_interaction[["Y_train"]] ; y_test = result_interaction[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = T;
  # s0 = result_interaction[["S0"]]; p = 1000 ;
  # model = "scad"; topgenes = NULL; true_beta = result_interaction[["beta_truth"]]
  
  # result[["clustersAddon"]] %>% print(nrows=Inf)
  # result[["clustersAddon"]][, table(cluster, module)]
  # result %>% names
  # stability = F; gene_groups = result[["clustersAll"]];
  # x_train = result[["X_train"]] ; x_test = result[["X_test"]];
  # y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  # filter = F; filter_var = F; include_E = T; include_interaction = T;
  # s0 = result[["S0"]]; p = p ;true_beta = result[["beta_truth"]]
  # model = "lasso"; summary = "pc"; topgenes = NULL; clust_type="clust"; nPC = 1
  
  # model: "scad", "mcp", "lasso", "elasticnet", "ridge"
  # filter: T or F based on univariate filter
  
  print(paste(model,"filter = ", filter, "filter_var = ",filter_var, "include_E = ", include_E, "include_interaction = ", include_interaction, sep = " "))
  
  if (include_E == F & include_interaction == T) stop("include_E needs to be
                                                      TRUE if you want to include
                                                      interactions")
  #   if (filter == F & include_interaction == T) stop("Interaction can only be run
  #                                                      if filter is TRUE.
  #                                                      This is to avoid exceedingly
  #                                                      large models")
  if (is.null(topgenes) & filter == T) stop("Argument topgenes is missing but
                                            filter is TRUE. You need to provide
                                            a filtered list of genes if filter
                                            is TRUE")
  
  #gene.names <- colnames(x_train)[which(colnames(x_train) %ni% "E")]
  
  # mars model
  mars_model <- switch(model,
                      MARS = earth::earth(x = as.matrix(x_train), 
                                          y = y_train, 
                                          keepxy = TRUE,
                                          pmethod = "forward",
                                          nk = 1000,
                                          degree = 2, 
                                          trace = 4))
                                          # ,
                                          # nfold = 10))
  
  # selected genes
  # coef(mars_model)
  # get.used.pred.names(mars_model)
  # 
  # plot(mars_model, which=1, col.rsq=0) # which=1 for Model Selection plot only (optional)
  # plot.earth.models(mars_model$cv.list, which=1)
  # plot(mars_model)
  # plot(mars_model, which=1,
  #      col.mean.infold.rsq="blue", col.infold.rsq="lightblue",
  #      col.grsq=0, col.rsq=0, col.vline=0, col.oof.vline=0)
  

  # ONLY Jaccard index can be calculated for MARS
  # since get.used.pred.names returns only the non-zero coefficients, 
  # we give a coefficient of 1 here so that the stability calculation works
  # because it takes non-zero coef estimates
  
  coefs <- data.frame(get.used.pred.names(mars_model), rep(1, length(get.used.pred.names(mars_model)))) %>%
    as.data.table(keep.rownames = FALSE) %>%
    magrittr::set_colnames(c("Gene","coef.est"))

  if (stability) {
    # remove intercept for stability measures
    return(coefs)
  } else {
    
    mars.S.hat <- get.used.pred.names(mars_model)
    mars.S.hat.interaction <- grep(":", mars.S.hat, value = T)
    mars.S.hat.main <- setdiff(mars.S.hat, mars.S.hat.interaction)
    
    mars.pred <- predict(mars_model, newdata = x_test, trace = 4)
    
    # True Positive Rate
    mars.TPR <- length(intersect(mars.S.hat, s0))/length(s0)
    
    # True negatives
    trueNegs <- setdiff(colnames(x_train), s0)
    
    # these are the terms which the model identified as zero
    modelIdentifyZero <- setdiff(colnames(x_train),mars.S.hat)
    
    # how many of the terms identified by the model as zero, were actually zero
    # use to calculate correct sparsity as defined by Witten et al in the 
    # Cluster Elastic Net paper Technometrics 2013
    C1 <- sum(modelIdentifyZero %in% trueNegs)
    C2 <- length(intersect(mars.S.hat, s0))
    correct_sparsity <- (C1 + C2)/(ncol(x_train))
    
    # this is from Interaction Screening for Ultrahigh Dimensional Data by ning hao and hao helen zhang
    true.interaction_names <- grep(":", s0, value = T)
    true.main_effect_names <- setdiff(s0, true.interaction_names)
    
    all.interaction_names <- grep(":", colnames(x_train), value = T)
    all.main_effect_names <- setdiff(colnames(x_train), all.interaction_names)
    
    true.negative_main_effects <- setdiff(all.main_effect_names, true.main_effect_names)
    true.negative_interaction_effects <- setdiff(all.interaction_names, true.interaction_names)
    
    (mars.correct_zeros_main_effects <- sum(setdiff(all.main_effect_names, mars.S.hat.main) %in% true.negative_main_effects)/ length(true.negative_main_effects))
    (mars.correct_zeros_interaction_effects <- sum(setdiff(all.interaction_names, mars.S.hat.interaction) %in% true.negative_interaction_effects)/ length(true.negative_interaction_effects))
    
    (mars.incorrect_zeros_main_effects <- sum(setdiff(all.main_effect_names, mars.S.hat) %in% true.main_effect_names)/ length(true.main_effect_names))
    (mars.incorrect_zeros_interaction_effects <- sum(setdiff(all.interaction_names, mars.S.hat.interaction) %in% true.interaction_names)/ length(true.interaction_names))    
    
    # False Positive Rate = FP/(FP + TN) = FP / True number of 0 coefficients
    (mars.FPR <- sum(mars.S.hat %ni% s0)/(sum(mars.S.hat %ni% s0) + sum(modelIdentifyZero %in% trueNegs)))
    
    # # False Positive Rate
    # mars.FPR <- sum(mars.S.hat %ni% s0)/(p - length(s0))    
    
    # Mean Squared Error
    (mars.mse <- crossprod(mars.pred - y_test)/length(y_test))
    
    # Root Mean Squared Error
    (mars.RMSE <- sqrt(crossprod(mars.pred - y_test)/length(y_test)))
    
    # Mean Squared Error Oracle
    # mars.mse.oracle <- crossprod(mars.pred.oracle - y_test)/length(y_test)
    
    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)
    
    # sqrt(mse_null)
    
    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    # mars.r2 <- (mse_null - mars.mse)/mse_null
    
    # mars.adj.r2 <- 1 - (1 - mars.r2)*(length(y_test) - 1)/(length(y_test) - length(mars.S.hat) - 1)
    
    # model error
    # identical(true_beta %>% rownames(),coefs[["Gene"]])
    # mars.model.error <- {(true_beta - coefs[["coef.est"]]) %>% t} %*% WGCNA::cor(x_test[,coefs[["Gene"]]]) %*% (true_beta - coefs[["coef.est"]])
    
    ls <- list(mars.mse = as.numeric(mars.mse),
               mars.RMSE = as.numeric(mars.RMSE),
               # mars.r2 = as.numeric(mars.r2),
               # mars.adj.r2 = as.numeric(mars.adj.r2), 
               mars.S.hat = length(mars.S.hat),
               mars.TPR = mars.TPR, 
               mars.FPR = mars.FPR, 
               # mars.relative.mse = mars.mse/mars.mse.oracle, 
               # mars.model.error = mars.model.error,
               correct_sparsity,
               mars.correct_zeros_main_effects,
               mars.correct_zeros_interaction_effects,
               mars.incorrect_zeros_main_effects,
               mars.incorrect_zeros_interaction_effects
    )
    names(ls) <- c(paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_RMSE"),
                   # paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_r2"),
                   # paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_adjr2"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectSparsity"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroMain"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_CorrectZeroInter"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroMain"),
                   paste0("mars_na_",model,ifelse(include_interaction,"_yes","_no"),"_IncorrectZeroInter"))
    return(ls)
  }
}


print("starting MARS with interaction")

mars_res <- mapply(mars_fun,
                  #model = c("scad","mcp","lasso","elasticnet"),
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

print("done MARS with interaction")






# Mars simulation ESL page 327 ---------------------------------------------------------

result[["X_train"]] %>% dim


# Scenario 1

X <- replicate(3, rnorm(100))
dimnames(X)[[2]] <- c("X1","X2","e")
Y <- pmax(X[,"X1"]-1,0) + pmax(X[,"X1"]-1,0) * pmax(X[,"X2"]-0.8,0) + 0.12 * X[,"e"]
DT <- cbind(Y,X)

fit1 <- earth::earth(x = X[,c("X1","X2")], 
                     y = Y,
                     keepxy = TRUE,
                     pmethod = "forward",
                     # nk = 1000,
                     degree = 2, 
                     trace = 4,
                     nfold = 10, ncross = 5)
summary(fit1)
plotmo(fit1, which = 3)
plot(fit1, which=1, col.rsq=0) # which=1 for Model Selection plot only (optional)
plot.earth.models(fit1$cv.list, which=1)
plot(fit1)
plot(fit1, which=1,
     col.mean.infold.rsq="blue", col.infold.rsq="lightblue",
     col.grsq=0, col.rsq=0, col.vline=0, col.oof.vline=0)


# Used in manuscript ------------------------------------------------------

hingeprod <- function(x1, x2) {
  0.1*(x1 + x2) + 4 * pmax(x1-0.01, 0) * pmax(x2-0.05, 0)
  # 10*sin(pi * x*y)
  # tan(pmax(x,0) * pmax(y,0))
  # 0.1*exp(4*x) + 4/(1+exp(-20*(y-0.5)))
}


linprod <- function(x1, x2) {
  0.1*(x1 + x2)
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

result[["S0"]]
x1 <- result[["X_train"]][, result[["S0"]][1:25]]
dim(x1)
u1 <- svd(x1)$u[,1]
# u1 <- apply(x1, 1, mean)

x2 <- result[["X_train"]][, result[["S0"]][26:50]]
dim(x2)
u2 <- svd(x2)$u[,1]
# u2 <- apply(x2, 1, mean)

reslin <- outseq(u1,u2, linprod)
reshinge <- outseq(u1,u2, hingeprod)
# persp(res[["y"]], res[["x"]], res[["z"]])
# persp(res[["y"]], res[["x"]], res[["z"]],
#       theta = 10, phi = 30, expand = 0.8, col = "lightblue", ltheta = 120, shade = .01, ticktype = "simple",    xlab = "X1", ylab = "X2", zlab = "(X1-1)_+ * (X2-0.8)_+")

savepdf("~/git_repositories/eclust-simulation-aug2016/sim3-persp.pdf")
par(mfrow=c(1,2))
drape.plot(reslin[["y"]], reslin[["x"]], reslin[["z"]], col=colorRampPalette(brewer.pal(9, "Reds"))(100),
           horizontal = FALSE, add.legend = F, zlim = range(reshinge[["z"]]), zlim2 = range(reshinge[["z"]]), 
           theta = 50, phi = 30, expand = 0.75, ltheta = 120,
           xlab = "U1", ylab = "U2", zlab = "Y")#, main = "E = 0")
drape.plot(reshinge[["y"]], reshinge[["x"]], reshinge[["z"]], col=colorRampPalette(brewer.pal(9, "Reds"))(100),
           horizontal = FALSE, add.legend = F, zlim = range(reshinge[["z"]]), zlim2 = range(reshinge[["z"]]), 
           theta = 50, expand = 0.75, phi = 30, ltheta = 120,
           xlab = "U1", ylab = "U2", zlab = "Y")#, main = "E = 1")
dev.off()








# use different range for color scale and persp plot
# setting of border omits the mesh lines

 drape.plot(res[["y"]], res[["x"]], res[["z"]], 
            col=terrain.colors(128), border=NA)

persp.withcol <- function(x,y,z,pal,nb.col,...,xlg=TRUE,ylg=TRUE) {
	colnames(z) <- y
	rownames(z) <- x

	nrz <- nrow(z)
	ncz <- ncol(z) 

	color <- pal(nb.col)
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	facetcol <- cut(zfacet, nb.col)
	par(xlog=xlg,ylog=ylg)
	persp(
		as.numeric(rownames(z)),
		as.numeric(colnames(z)),
		as.matrix(z),
		col=color[facetcol],
		...)
}


library(fields)
data( volcano)
M<- nrow( volcano)
N<- ncol( volcano)
x<- seq( 0,1,,M)
y<- seq( 0,1,,N)

drape.plot( x,y,volcano, col=terrain.colors(128))-> pm 

# use different range for color scale and persp plot
# setting of border omits the mesh lines

 drape.plot( x,y,volcano, col=terrain.colors(128),zlim=c(0,300),
                     zlim2=c( 120,165), border=NA)















library(plot3D)
surf3D(res[["z"]], res[["x"]], res[["y"]], colvar = z, colkey = FALSE, facets = FALSE)

persp3D(z = res[["z"]], phi = 30,
colkey = list(length = 0.2, width = 0.4, shift = 0.15,
cex.axis = 0.8, cex.clab = 0.85), lighting = TRUE, lphi = 90,
clab = c("","height","m"), bty = "f")
volcano

persp(res[["x"]], res[["y"]], res[["z"]], theta = 135, phi = 30, col = "green3", scale = FALSE,
      ltheta = -120, shade = 0.75)

x <- seq(-10, 10, length= 30)
y <- x
f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
z <- outer(x, y, f)
z[is.na(z)] <- 1
op <- par(bg = "white")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      ltheta = 120, shade = 0.75, ticktype = "detailed",
      xlab = "X", ylab = "Y", zlab = "Sinc( r )"
) -> res
round(res, 3)

