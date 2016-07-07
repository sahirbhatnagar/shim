##################################
# R source code file for analyzing NIHPD data
# Created by Sahir, June 20
# Updated:
# NOTE: 
# 
# 
##################################

rm(list=ls())
dev.off()
source("functions.R")
source("packages.R")
source("data_cleaning.R")
Rcpp::sourceCpp("kmeans.cpp")
options(scipen = 999, digits = 4)


# get clusters of data to reduce the dimension of the prooblem ------------

phenotypeVariable <- "WASI_Full_Scale_IQ" 
# phenotypeVariable <- "WASI_Verbal_IQ"
# phenotypeVariable <- "WASI_Performance_IQ"
# exposureVariable <- "E"
exposureVariable <- "income_binary2"

# this is the cortical thickness data for the 1st timepoint only
# all the data for all timepoints is in dat$Thick.81924
thicknessMat <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, brain_probes, with = F] %>% as.matrix()

# DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, brain_probes, with = F] %>% as.matrix() %>% dim

kclustData <- mlKmeans(thicknessMat, 50)

table(kclustData$result)

# these contain the brain probes in their respective kmeans clusters
# subjects are rows, columns are probes
clusters <- lapply(by(t(thicknessMat),kclustData$result,identity), function(i) t(as.matrix(i)))
lapply(clusters, dim)
lapply(clusters, is.matrix)

# lapply(clusters, function(i) all(colnames(i %in% brain_probes))) %>% unlist() %>% all


# Create training and test set indices ------------------------------------

Y <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, phenotypeVariable, with = F] %>% as.matrix() %>% drop

trainIndex <- caret::createDataPartition(Y, p = 2/3, list = FALSE, times = 1) %>% drop
testIndex <- which(seq_len(length(Y)) %ni% trainIndex)
# all(c(trainIndex, testIndex) %in% seq_len(length(Y)))

# we also need the exposure variable as a vector
E <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, exposureVariable, with = F] %>% as.matrix() %>% drop




# PLS ---------------------------------------------------------------------

# library(pls)
# 
# 
# xy_train <- data.frame(Y = as.matrix(Y[trainIndex]), clusters[[48]][trainIndex,])
# xy_test <- data.frame(Y = as.matrix(Y[testIndex]), clusters[[48]][testIndex,])
# 
# gas1 <- pls::plsr(Y ~ ., ncomp = 2, data = xy_train, validation = "CV")
# summary(gas1)
# explvar(gas1)
# plot(RMSEP(gas1), legendpos = "topright")
# plot(gas1, ncomp = 1, asp = 1, line = TRUE)
# plot(gas1, "loadings", comps = 1, legendpos = "topleft",
#      xlab = "nm")
# abline(h = 0)
# predict(gas1, ncomp = 2, type = "scores")
# 
# predict(gas1, ncomp = 1, newdata = xy_test, type = "scores")[,1]

# This is where the loop begins -------------------------------------------
# loop over each element in 'clusters'

cluster_kmeans <- function(data, 
                        response, 
                        exposure, 
                        train_index, 
                        test_index,
                        min_cluster_size = 50,
                        cluster_distance = c("corr", "corr0", "corr1", "tom",
                                             "tom0", "tom1", "diffcorr",
                                             "difftom", "fisherScore"),
                        cluster_method = c("hclust", "protoclust"),
                        cut_method = c("dynamic", "gap", "fixed"),
                        agglomeration_method = c("average", "complete", "ward.D2",
                                   "single", "ward.D", "mcquitty",
                                   "median", "centroid"),
                        distance_method = c("euclidean", "maximum", "manhattan",
                                           "canberra", "binary", "minkowski"),
                        eclust_add = TRUE, 
                        eclust_distance = c("fisherScore", "corScor", "diffcorr",
                                              "difftom"),
                        nPC = 1,
                        n_clusters, K.max, B) {
  
  # data <- clusters[[1]]  
  # exposure <- E
  # train_index <- trainIndex
  # test_index <- testIndex
  # numberPC = 1
  # cluster_distance <- "tom"
#   cluster_method <- "hclust"
#   cut_method <- "dynamic"
#   agglomeration_method <- "average"
#   distance_method <- "euclidean"
#   eclust_add <- TRUE 
#   eclust_distance <- "difftom"
#   cut_method <- "dynamic" #,"gap", "fixed")
  # summary <- match.arg(summary)
  
  # ==================================================================
  
  xtrain <- data[train_index,]
  xtest <- data[test_index,]
  etrain <- exposure[train_index]
  
  ytrain <- response[train_index]
  ytest <- response[test_index]
  
  xtrainE0 <- xtrain[which(etrain == 0), ]
  xtrainE1 <- xtrain[which(etrain == 1), ]

  corrTrain <- WGCNA::cor(xtrain)

  message(sprintf("Calculating number of clusters based on %s using %s with %s
  linkage and the %s method to determine the number of clusters",
                  cluster_distance, cluster_method, agglomeration_method, cut_method))
  
  #############################################################################
  #               CORRELATION/TOM CLUSTERS                                    #
  #############################################################################
  
  # clusters based on cluster_distance argument
  similarity <- switch(cluster_distance,
                       corr = corrTrain,
                       corr0 = WGCNA::cor(xtrainE0),
                       corr1 = WGCNA::cor(xtrainE1),
                       diffcorr = {
                         corr0 <- WGCNA::cor(xtrainE0)
                         corr1 <- WGCNA::cor(xtrainE1)
                         abs(corr1 - corr0)
                         },
                       difftom = {
                         tomTrainE0 <- WGCNA::TOMsimilarityFromExpr(xtrainE0)
                         tomTrainE1 <- WGCNA::TOMsimilarityFromExpr(xtrainE1)
                         tomTrainDiff <- abs(tomTrainE1 - tomTrainE0)
                         dimnames(tomTrainDiff)[[1]] <- dimnames(corrTrain)[[1]]
                         dimnames(tomTrainDiff)[[2]] <- dimnames(corrTrain)[[2]]
                         tomTrainDiff
                       },
                       tom0 = {
                         tomTrainE0 <- WGCNA::TOMsimilarityFromExpr(xtrainE0)
                         dimnames(tomTrainE0)[[1]] <- dimnames(corrTrain)[[1]]
                         dimnames(tomTrainE0)[[2]] <- dimnames(corrTrain)[[2]]
                         tomTrainE0
                       },
                       tom1 = {
                         tomTrainE1 <- WGCNA::TOMsimilarityFromExpr(xtrainE1)
                         dimnames(tomTrainE1)[[1]] <- dimnames(corrTrain)[[1]]
                         dimnames(tomTrainE1)[[2]] <- dimnames(corrTrain)[[2]]
                         tomTrainE1
                       },
                       tom = {
                         tomTrainAll <- WGCNA::TOMsimilarityFromExpr(xtrain)
                         dimnames(tomTrainAll)[[1]] <- dimnames(corrTrain)[[1]]
                         dimnames(tomTrainAll)[[2]] <- dimnames(corrTrain)[[2]]
                         tomTrainAll
                       },
                       fisherScore = {
                         n0 <- nrow(xtrainE0)
                         n1 <- nrow(xtrainE1)
                         corr0 <- WGCNA::cor(xtrainE0)
                         corr1 <- WGCNA::cor(xtrainE1)
                         fisherZ(n0 = n0, cor0 = corr0,
                                 n1 = n1, cor1 = corr1)
                       })
  
  # results for clustering
  # note that the cluster_similarity returns the PCs but you dont need this info
  # because it is calculated in the clust_fun fitting function
  # we just need to provide the clust_fun function the group membership for all the data
  res <- if (cluster_distance %in% c("diffcorr", "difftom", "fisherScore")) {
    cluster_similarity(x = similarity,
                       x_train = xtrain,
                       x_test = xtest,
                       y_train = ytrain,
                       y_test = ytest,
                       distance_method = distance_method,
                       cluster_method = cluster_method,
                       cut_method = cut_method,
                       agglomeration_method = agglomeration_method,
                       K.max = K.max, B = B, 
                       n_clusters = n_clusters,
                       min_cluster_size = min_cluster_size,
                       nPC = nPC)
  } else {
    cluster_similarity(x = similarity,
                       x_train = xtrain,
                       x_test = xtest,
                       y_train = ytrain,
                       y_test = ytest,
                       min_cluster_size = min_cluster_size,
                       cluster_method = cluster_method,
                       cut_method = cut_method,
                       agglomeration_method = agglomeration_method,
                       K.max = K.max, B = B, n_clusters = n_clusters, nPC = nPC)
  }
  
  #############################################################################
  #               ECLUST CLUSTERS                                             #
  #############################################################################
  
  message(paste("Calculating number of environment clusters based on",eclust_distance))
  
  # clusters based on eclust_distance
  similarityEclust <- switch(eclust_distance,
                             corr = corrTrain,
                             corr0 = WGCNA::cor(xtrainE0),
                             corr1 = WGCNA::cor(xtrainE1),
                             diffcorr = {
                               corr0 <- WGCNA::cor(xtrainE0)
                               corr1 <- WGCNA::cor(xtrainE1)
                               abs(corr1 - corr0)
                             },
                             difftom = {
                               tomTrainE0 <- WGCNA::TOMsimilarityFromExpr(xtrainE0)
                               tomTrainE1 <- WGCNA::TOMsimilarityFromExpr(xtrainE1)
                               tomTrainDiff <- abs(tomTrainE1 - tomTrainE0)
                               dimnames(tomTrainDiff)[[1]] <- dimnames(corrTrain)[[1]]
                               dimnames(tomTrainDiff)[[2]] <- dimnames(corrTrain)[[2]]
                               tomTrainDiff
                             },
                             tom0 = {
                               tomTrainE0 <- WGCNA::TOMsimilarityFromExpr(xtrainE0)
                               dimnames(tomTrainE0)[[1]] <- dimnames(corrTrain)[[1]]
                               dimnames(tomTrainE0)[[2]] <- dimnames(corrTrain)[[2]]
                               tomTrainE0
                             },
                             tom1 = {
                               tomTrainE1 <- WGCNA::TOMsimilarityFromExpr(xtrainE1)
                               dimnames(tomTrainE1)[[1]] <- dimnames(corrTrain)[[1]]
                               dimnames(tomTrainE1)[[2]] <- dimnames(corrTrain)[[2]]
                               tomTrainE1
                             },
                             tom = {
                               tomTrainAll <- WGCNA::TOMsimilarityFromExpr(xtrain)
                               dimnames(tomTrainAll)[[1]] <- dimnames(corrTrain)[[1]]
                               dimnames(tomTrainAll)[[2]] <- dimnames(corrTrain)[[2]]
                               tomTrainAll
                             },
                             fisherScore = {
                               n0 <- nrow(xtrainE0)
                               n1 <- nrow(xtrainE1)
                               corr0 <- WGCNA::cor(xtrainE0)
                               corr1 <- WGCNA::cor(xtrainE1)
                               fisherZ(n0 = n0, cor0 = corr0,
                                       n1 = n1, cor1 = corr1)
                             })
  
  resEclust <- if (eclust_distance %in% c("diffcorr","difftom",
                                            "corScor", "tomScor","fisherScore")) {
    cluster_similarity(x = similarityEclust,
                       x_train = xtrain,
                       x_test = xtest,
                       y_train = ytrain,
                       y_test = ytest,
                       min_cluster_size = min_cluster_size,
                       distance_method = distance_method,
                       cluster_method = cluster_method,
                       cut_method = cut_method,
                       agglomeration_method = agglomeration_method,
                       K.max = K.max, B = B, n_clusters = n_clusters, nPC = nPC)
  } else {
    cluster_similarity(x = similarityEclust,
                       x_train = xtrain,
                       x_test = xtest,
                       y_train = ytrain,
                       y_test = ytest,
                       min_cluster_size = min_cluster_size,
                       cluster_method = cluster_method,
                       cut_method = cut_method,
                       agglomeration_method = agglomeration_method,
                       K.max = K.max, B = B, n_clusters = n_clusters, nPC = nPC)
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
  
  # this contains the clusters from the cluster_distance (e.g. TOM matrix)
  # and the clusters from the eclust_distance (e.g. diffTOM)
  clustersAddon <- rbindlist(list(clustersAll, clustersEclust))
  # clustersAddon[, table(cluster, module)]

  # these are only derived on the main effects genes.. E is only included in the model
  # this is the clusters based on tom and difftom
  PC_and_avg_Addon <- extractPC(x_train = xtrain[,clustersAddon$gene],
                                colors = clustersAddon$cluster,
                                x_test = xtest[,clustersAddon$gene],
                                y_train = ytrain,
                                y_test = ytest,
                                nPC = nPC)
  
  # this is the clusters based on tom only
  PC_and_avg_All <- extractPC(x_train = xtrain[,clustersAll$gene],
                              colors = clustersAll$cluster,
                              x_test = xtest[,clustersAll$gene],
                              y_train = ytrain,
                              y_test = ytest,
                              nPC = nPC)

  # n.clusters <- PC_and_avg_Addon$nclusters
  
  # this contains either the averages or PCs for each module in a data.frame
#   clust_data_Addon <- switch(summary,
#                        avg = PC_and_avg_Addon$averageExpr,
#                        pc = PC_and_avg_Addon$PC)
#   
#   clust_data_All <- switch(summary,
#                            avg = PC_and_avg_All$averageExpr,
#                            pc = PC_and_avg_All$PC)
  
  
  return(list(clustersAddon = PC_and_avg_Addon, 
              clustersAll = PC_and_avg_All, 
              etrain = etrain))
}


# pl <- lapply(clusters[c(4,5)], function(i) cluster_kmeans(data = i, 
#                                                     exposure = E,
#                                                     response = Y,
#                                                     min_cluster_size = 50,
#                                                     train_index = trainIndex, 
#                                                     test_index = testIndex,
#                                                     cluster_distance = "tom",
#                                                     cluster_method = "hclust",
#                                                     cut_method = "dynamic",
#                                                     agglomeration_method = "average",
#                                                     distance_method = "euclidean",
#                                                     eclust_add = TRUE,
#                                                     eclust_distance = "difftom",
#                                                     nPC = 1))

doMC::registerDoMC(cores = 30)
pp <- mclapply(clusters, function(i) cluster_kmeans(data = i, 
                                                 exposure = E,
                                                 response = Y,
                                                 min_cluster_size = 50,
                                                 train_index = trainIndex, 
                                                 test_index = testIndex,
                                                 cluster_distance = "tom",
                                                 cluster_method = "hclust",
                                                 cut_method = "dynamic",
                                                 agglomeration_method = "average",
                                                 distance_method = "euclidean",
                                                 eclust_add = TRUE,
                                                 eclust_distance = "difftom",
                                                 nPC = 1),
               mc.cores = 30)

save(pp, file = "PC_NIHPD_based_on_all_data_converted_to_50_kmeans_clusters.RData")

# pp$`0`$clustersAddon$PLS
# pp$`44`$clustersAddon$varExplained %>% plot
# pp$`43`$clustersAddon$varExplained %>% plot
# pheatmap::pheatmap(pp$`44`$clustersAddon$PC[which(pp$`44`$etrain==0),])
# pheatmap::pheatmap(pp$`44`$clustersAddon$PC[which(pp$`44`$etrain==1),])

# to combine all the principal components
load("PC_NIHPD_based_on_all_data_converted_to_50_kmeans_clusters.RData")
pcTrain <- do.call(cbind, lapply(pp, function(i) i[["clustersAddon"]][["PC"]]))
avgTrain <- do.call(cbind, lapply(pp, function(i) i[["clustersAddon"]][["averageExpr"]]))
colnames(pcTrain) <- gsub("\\.","_",colnames(pcTrain))
colnames(pcTrain) <- paste0("PC",colnames(pcTrain))
datt <- cbind(pcTrain, Y = Y[trainIndex])
colnames(datt)
str(datt)

dim(pcTrain)
pcTest <- do.call(cbind, lapply(pp, function(i) i[["clustersAddon"]][["PCTest"]]))
avgTest <- do.call(cbind, lapply(pp, function(i) i[["clustersAddon"]][["averageExprTest"]]))
dim(pcTest)


# to get number of clusters from each Kmeans cluster
do.call(c, lapply(pp, function(i) i[["clustersAddon"]][["nclusters"]])) %>% plot

# variance explained
do.call(c, lapply(pp, function(i) i[["clustersAddon"]][["varExplained"]])) %>% plot

library(ComplexHeatmap)
require(circlize)
pcTrain %>% dim
max(pcTrain);min(pcTrain)
apply(pcTrain, 2, max) %>% as.numeric() %>% max
apply(pcTrain, 2, min) %>% as.numeric() %>% min


cm <- colorRamp2(seq(min(pcTrain), max(pcTrain), length.out = 50), viridis(50))
ht1 = Heatmap(t(pcTrain[which(pp[[1]][["etrain"]]==0),]), 
              name = "E=0",
              # col = viridis(100),
              col = cm,
              # column_title = "E = 0 : Age [4.8, 11.3]",
              column_title = "Income_Level: 1-7",
              show_row_names = FALSE)
ht2 = Heatmap(t(pcTrain[which(pp[[1]][["etrain"]]==1),]), 
              name = "E=1",
              # col = viridis(100), 
              col = cm,
              # column_title = "E = 1 : Age [11.3, 18]",
              column_title = "Income_Level: 8-10",
              show_row_names = FALSE)
ht1 + ht2


cv.fit <- cv.glmnet(x = as.matrix(pcTrain), y = Y[trainIndex], alpha = 0.5, standardize = T, intercept=T)
cv.fit <- cv.glmnet(x = as.matrix(avgTrain), y = Y[trainIndex], alpha = 0.5)

plot(cv.fit)
as.matrix(coef(cv.fit, s = "lambda.1se"))[which(as.matrix(coef(cv.fit, s = "lambda.1se"))!=0),,drop=F]
as.matrix(coef(cv.fit, s = "lambda.min"))[which(as.matrix(coef(cv.fit, s = "lambda.min"))!=0),,drop=F]
cv.fit$cvm[which.min(cv.fit$cvm)]




res <- ldply(colnames(pcTrain), function(i) {
  
  # i = "49.pc1_4"
  # paste(sprintf("`%s`", i), collapse="+")
  include_E = FALSE
  include_interaction = FALSE
  fit <- lm(as.formula(paste("Y ~",sprintf("`%s`", i), if (include_E) "+E", if (include_interaction) paste0("+",sprintf("`%s`", i),":E"))), data = datt)
  #fit %>% summary
  # fit$model
  coef.index <- if (include_interaction) paste0(sprintf("`%s`", i),":E") else sprintf("`%s`", i)
  
  data.frame(name = i, pvalue = as.numeric(pt(abs(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
                                    df = fit$df.residual, lower.tail = F)*2),
             test.stat = as.numeric(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
             'mean' = mean(datt[,i]),
             'sd' = sd(datt[,i]))
})

res[order(res$pvalue)[1:10],]
hist(res$pvalue)

manhattanly::qqly(res, p = "pvalue")

qvalue::qvalue(res$pvalue)$qvalues[order(qvalue::qvalue(res$pvalue)$qvalues)[1:100]]



clusters[[31]] %>% dim






cv.fit2 <- cv.glmnet(x = thicknessMat[trainIndex,], y = Y[trainIndex], alpha = 0.5)
plot(cv.fit2)
as.matrix(coef(cv.fit2, s = "lambda.1se"))[which(as.matrix(coef(cv.fit2, s = "lambda.1se"))!=0),,drop=F]
as.matrix(coef(cv.fit2, s = "lambda.min"))[which(as.matrix(coef(cv.fit2, s = "lambda.min"))!=0),,drop=F]
cv.fit2$cvm[which.min(cv.fit2$cvm)]



crossprod(as.matrix(Y[testIndex]) - predict(cv.fit, newx = as.matrix(pcTest), s = "lambda.min"))
crossprod(as.matrix(Y[testIndex]) - predict(cv.fit2, newx = thicknessMat[testIndex,], s = "lambda.min"))



# Create Interaction Data for cv.shim ----------------------------------------

prepare_data <- function(data, response = "Y", exposure = "E", probe_names) {
  
  # data = cbind(pcTrain, Y = Y[trainIndex], E = E[trainIndex])
  
    
  # ===========================================================
    
  # Check for sensible dataset
  ## Make sure you have response, exposure.
  if (!(response %in% colnames(data))) stop(sprintf("response argument specified as %s but this column not found in 'data' data.frame", response))
  if (!(exposure %in% colnames(data))) stop(sprintf("exposure argument specified as %s but this column not found in 'data' data.frame", exposure))
  if (!missing(probe_names)) {
    if (!(probe_names %in% colnames(data))) stop(sprintf("probe_names argument specified as %s but this column not found in 'data' data.frame", probe_names))
  }
  
  # if missing main_effect_names, assume everything except response and exposure 
  # are main effects
  if (missing(probe_names)) {
    probe_names <- setdiff(colnames(data), c(response, exposure))
  }
  
  # rename response to be Y and exposure to be E
  colnames(data)[which(colnames(data) == response)] <- "Y"
  colnames(data)[which(colnames(data) == exposure)] <- "E"
  
  x_mat <- model.matrix(as.formula(paste0("~(", paste0(probe_names, collapse="+"), ")*E - 1")), data = data)
  
  
  # reformulate(paste0("~(", paste0(colnames(pcTrain)[1:5], collapse="+"), ")*E"))
  
  
  interaction_names <- grep(":", colnames(x_mat), value = T)
  main_effect_names <- setdiff(colnames(x_mat), interaction_names)
  
  return(list(X = x_mat, Y = data[["Y"]], E = data[["E"]], 
              main_effect_names = main_effect_names, 
              interaction_names = interaction_names))
  # x_mat
  
  
}

kl <- prepare_data(data = cbind(pcTrain, pheno = Y[trainIndex], income = E[trainIndex]),
                   response = "pheno", exposure = "income")

kl$X %>% dim()
kl$main_effect_names
kl$interaction_names

cv.fit2 <- cv.glmnet(x = kl$X, y = kl$Y, alpha = 0.5)
plot(cv.fit2)
as.matrix(coef(cv.fit2, s = "lambda.1se"))[which(as.matrix(coef(cv.fit2, s = "lambda.1se"))!=0),,drop=F]
as.matrix(coef(cv.fit2, s = "lambda.min"))[which(as.matrix(coef(cv.fit2, s = "lambda.min"))!=0),,drop=F]

library(eclust)
library(doMC)
registerDoMC(cores = 10)
system.time(cv_shim <- cv.shim(x =kl$X, y = kl$Y,
                   main.effect.names = kl$main_effect_names,
                   interaction.names = kl$interaction_names,
                   parallel = TRUE, verbose = TRUE,
                   type.measure = c("mse"), 
                   nfolds = 10))

library(ggplot2)
plot(cv_shim)
as.matrix(coef(cv_shim$shim.fit, s = "lambda.1se"))[which(as.matrix(coef(cv_fit$shim.fit, s = "lambda.1se"))!=0),,drop=F]


coef(cv_shim$shim.fit, s = "lambda.1se")



