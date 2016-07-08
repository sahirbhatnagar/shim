##################################
# R source code file for analyzing bouchard data on qsub
# Created by Sahir, July 4, 2016
# Updated:
# NOTE: Hosted on 'bouchard' branch of eclust repo on github.
# This file is similar to the clust_kmeans_bouchard.R file on the bouchard branch of 
# the eclust repo
# this file is also the one being used for the analysis.
# The 'analysis.R' was just used previously for implementing cor_scor
# 
##################################


source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/packages.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/data_cleaning.R")

Rcpp::sourceCpp("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/kmeans.cpp")
options(scipen = 999, digits = 4)

load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/kclustData.RData")

# these contain the cpg probes in their respective kmeans clusters
# subjects are rows, columns are probes
clusters <- lapply(by(placentaALL,kclustData$result,identity), function(i) t(as.matrix(i)))
lapply(clusters, dim)

parameterIndex <- as.numeric(commandArgs(trailingOnly = T))
parameterIndex

# Create training and test set indices ------------------------------------

Y <- DT.pheno.placenta$`imc z-score` %>% as.matrix() %>% drop

trainIndex <- caret::createDataPartition(Y, p = 1, list = FALSE, times = 1) %>% drop

# not really used... need to do LOOCV 
# testIndex <- which(seq_len(length(Y)) %ni% trainIndex)
testIndex <- trainIndex
# all(c(trainIndex, testIndex) %in% seq_len(length(Y)))

# we also need the exposure variable as a vector
# columns of methylation data match the ordering of the phenotype data
# see 'data_cleaning.R' for details
# needs to be 0 and 1
E <- as.numeric(DT.pheno.placenta$case)-1

# This is where the loop begins -------------------------------------------
# loop over each element in 'clusters'

doMC::registerDoMC(cores = 20)
pp <- mclapply(clusters[parameterIndex], function(i) cluster_kmeans(data = i, 
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
               mc.cores = 20)



save(pp, file = paste0(
  "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git",tempfile(tmpdir = ""),"PC_bouchard_based_on_all_data_converted_to_50_kmeans_clusters.RData"
  )
)


rm(pl)
load("file75862e2e3908PC_bouchard_based_on_all_data_converted_to_50_kmeans_clusters.RData")
c1 <- pp
rm(pp)
load("file336f78f6664cPC_bouchard_based_on_all_data_converted_to_50_kmeans_clusters.RData")
c2 <- pp
rm(pp)
pp <- c(c1,c2)


pp[[1]] %>% str


library(rlist)
pp_filter <- list.filter(pp, clustersAddon$nclusters < 500)


# to get number of clusters from each Kmeans cluster
N_clusters <- do.call(rbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["nclusters"]]))

do.call(c, lapply(pp_filter, function(i) i[["clustersAddon"]][["nclusters"]])) %>% plot

# variance explained
do.call(c, lapply(pp_filter, function(i) i[["clustersAddon"]][["varExplained"]])) %>% plot

# to combine all the principal components
pcTrain <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["PC"]]))
avgTrain <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["averageExpr"]]))
colnames(pcTrain) <- gsub("\\.","_",colnames(pcTrain))
colnames(pcTrain) <- paste0("PC",colnames(pcTrain))
datt <- cbind(pcTrain, Y = Y[trainIndex])
colnames(datt)
str(datt)

dim(pcTrain)
pcTest <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["PCTest"]]))
avgTest <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["averageExprTest"]]))
dim(pcTest)

library(ComplexHeatmap)
require(circlize)

cm <- colorRamp2(seq(min(pcTrain), max(pcTrain), length.out = 100), viridis(100))
ht1 = Heatmap(t(pcTrain[which(pp_filter[[1]][["etrain"]]==0),]), 
              name = "E=0",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 0 : Age [4.8, 11.3]",
              # column_title = "Income_Level: 1-7",
              column_title = "NGD",
              show_row_names = FALSE)
ht2 = Heatmap(t(pcTrain[which(pp_filter[[1]][["etrain"]]==1),]), 
              name = "E=1",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 1 : Age [11.3, 18]",
              # column_title = "Income_Level: 8-10",
              column_title = "GD",
              show_row_names = FALSE)
ht1 + ht2

cm <- colorRamp2(seq(min(avgTrain), max(avgTrain), length.out = 100), viridis(100))
ht1 = Heatmap(t(avgTrain[which(pp_filter[[1]][["etrain"]]==0),]), 
              name = "E=0",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 0 : Age [4.8, 11.3]",
              # column_title = "Income_Level: 1-7",
              column_title = "NGD",
              show_row_names = FALSE)
ht2 = Heatmap(t(avgTrain[which(pp_filter[[1]][["etrain"]]==1),]), 
              name = "E=1",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 1 : Age [11.3, 18]",
              # column_title = "Income_Level: 8-10",
              column_title = "GD",
              show_row_names = FALSE)
ht1 + ht2








cv.fit <- cv.glmnet(x = as.matrix(pcTrain), y = Y[trainIndex], alpha = 0.5, standardize = T, intercept=T)
# cv.fit <- cv.glmnet(x = as.matrix(avgTrain), y = Y[trainIndex], alpha = 0.5)

plot(cv.fit)
as.matrix(coef(cv.fit, s = "lambda.1se"))[which(as.matrix(coef(cv.fit, s = "lambda.1se"))!=0),,drop=F]
as.matrix(coef(cv.fit, s = "lambda.min"))[which(as.matrix(coef(cv.fit, s = "lambda.min"))!=0),,drop=F]
cv.fit$cvm[which.min(cv.fit$cvm)]


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

# cv.fit2 <- cv.glmnet(x = kl$X, y = kl$Y, alpha = 0.5)
# plot(cv.fit2)
# as.matrix(coef(cv.fit2, s = "lambda.1se"))[which(as.matrix(coef(cv.fit2, s = "lambda.1se"))!=0),,drop=F]
# as.matrix(coef(cv.fit2, s = "lambda.min"))[which(as.matrix(coef(cv.fit2, s = "lambda.min"))!=0),,drop=F]

library(eclust)
library(doMC)
registerDoMC(cores = 10)
system.time(shim <- shim(x =kl$X, y = kl$Y,
                               main.effect.names = kl$main_effect_names,
                               interaction.names = kl$interaction_names,
                               verbose = TRUE))

system.time(cv_shim <- cv.shim(x =kl$X, y = kl$Y,
                               main.effect.names = kl$main_effect_names,
                               interaction.names = kl$interaction_names,
                               parallel = TRUE, verbose = TRUE,
                               type.measure = c("mse"), 
                               nfolds = 10))

