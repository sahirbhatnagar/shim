##################################
# R source code file for analyzing NIHPD data using qsub
# Created by Sahir, July 7
# Updated:
# NOTE: 
# 
# 
##################################

# rm(list=ls())
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_git/functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_git/packages.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_git/data_cleaning.R")
Rcpp::sourceCpp("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_git/kmeans.cpp")
options(scipen = 999, digits = 4)



# Filter out probes -------------------------------------------------------

# Response is WASI_Full_Scale_IQ,  p-values per probe based on 
# lm controlling for Subject_Gender + Site_Location
load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_git/uni_results.RData")

# uni_results[, cv:=sd/mean, by = probe]
# 
# theta = seq(from=0, to=0.99, by=0.02)
# filterChoices = data.frame(
#   `mean` = uni_results$mean,
#   `geneID` = 1:nrow(uni_results),
#   `sd` = uni_results$sd,
#   `cv` = uni_results$cv
# )
# 
# rejChoices = sapply(filterChoices, function(f)
#   filtered_R(alpha=0.1, filter=f, test=uni_results$pvalue, theta=theta, method="BH"))

# library("RColorBrewer")
# myColours = brewer.pal(ncol(filterChoices), "Set1")
# matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
#         xlab=expression(theta), ylab="number of rejections")
# legend("topleft", legend=colnames(filterChoices), fill=myColours)
# abline(v = theta[order(rejChoices[,"mean"], decreasing = TRUE)][1], lty=2)
# 
# theta[order(rejChoices[,"mean"], decreasing = TRUE)]
# quantile(uni_results$mean)

# thetaThreshold <- theta[order(rejChoices[,"mean"], decreasing = TRUE)][1]
# filtered_probes <- uni_results[with(uni_results, mean > quantile(mean, thetaThreshold)),]$probe

# thetaThreshold <- theta[order(rejChoices[,"mean"], decreasing = TRUE)][1]
filtered_probes <- uni_results[with(uni_results, sd > quantile(sd, 0.877946)),]$probe
# plot(uni_results[probe %in% filtered_probes]$sd,ylim=c(0,1))
# filter based on sd

# head(filtered_probes)


# get clusters of data to reduce the dimension of the prooblem ------------

phenotypeVariable <- "WASI_Full_Scale_IQ" 
# phenotypeVariable <- "WASI_Verbal_IQ"
# phenotypeVariable <- "WASI_Performance_IQ"
# exposureVariable <- "E"
exposureVariable <- "income_binary2"

# this is the cortical thickness data for the 1st timepoint only
# all the data for all timepoints is in dat$Thick.81924
thicknessMat <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, brain_probes, with = F] %>% as.matrix()
thicknessMat <- thicknessMat[,filtered_probes]
thicknessMat %>% dim

# DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, brain_probes, with = F] %>% as.matrix() %>% dim

# kclustData <- mlKmeans(thicknessMat, 50)
# save(kclustData, file = "kclustData_NIH_50_kmeans_after_mean_filter.RData")

# load("kclustData_NIH_50_kmeans_after_mean_filter.RData")
# table(kclustData$result)

# these contain the brain probes in their respective kmeans clusters
# subjects are rows, columns are probes
# clusters <- lapply(by(t(thicknessMat),kclustData$result,identity), function(i) t(as.matrix(i)))
# lapply(clusters, dim)
# lapply(clusters, is.matrix)

# lapply(clusters, function(i) all(colnames(i %in% brain_probes))) %>% unlist() %>% all

# Create training and test set indices ------------------------------------

Y <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, phenotypeVariable, with = F] %>% as.matrix() %>% drop

trainIndex <- caret::createDataPartition(Y, p = 1, list = FALSE, times = 1) %>% drop
testIndex <- trainIndex
# testIndex <- which(seq_len(length(Y)) %ni% trainIndex)
# all(c(trainIndex, testIndex) %in% seq_len(length(Y)))

# we also need the exposure variable as a vector
E <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, exposureVariable, with = F] %>% as.matrix() %>% drop


# doMC::registerDoMC(cores = 5)
# pp <- mclapply(clusters, function(i) cluster_kmeans(data = i, 
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
#                                                     nPC = 1),
#                mc.cores = 5)



pp <- cluster_kmeans(data = thicknessMat, 
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
               nPC = 1)

save(pp, file = "PC_NIHPD_10000_most_variable_TOM_DIFFTOM.RData")


# doMC::registerDoMC(cores = 5)
# pp2 <- mclapply(clusters, function(i) cluster_kmeans(data = i, 
#                                                     exposure = E,
#                                                     response = Y,
#                                                     min_cluster_size = 50,
#                                                     train_index = trainIndex, 
#                                                     test_index = testIndex,
#                                                     cluster_distance = "corr",
#                                                     cluster_method = "hclust",
#                                                     cut_method = "dynamic",
#                                                     agglomeration_method = "average",
#                                                     distance_method = "euclidean",
#                                                     eclust_add = TRUE,
#                                                     eclust_distance = "fisherScore",
#                                                     nPC = 1),
#                mc.cores = 5)

pp2 <- cluster_kmeans(data = thicknessMat, 
                      exposure = E,
                      response = Y,
                      min_cluster_size = 50,
                      train_index = trainIndex, 
                      test_index = testIndex,
                      cluster_distance = "corr",
                      cluster_method = "hclust",
                      cut_method = "dynamic",
                      agglomeration_method = "average",
                      distance_method = "euclidean",
                      eclust_add = TRUE,
                      eclust_distance = "fisherScore",
                      nPC = 1)

# save(pp2, file = "PC_NIHPD_based_on_all_data_converted_to_50_kmeans_clusters_filtered_mean_corr_fisherScore.RData")
save(pp2, file = "PC_NIHPD_10000_most_variable_corr_fisherScore.RData")





# after submitting to qsub ------------------------------------------------

# rm(pp_filter,pp, pp2)
# load("PC_NIHPD_based_on_all_data_converted_to_50_kmeans_clusters_filtered_mean_TOM_DIFFTOM.RData")
# 
# # load("PC_NIHPD_based_on_all_data_converted_to_50_kmeans_clusters_filtered_mean_corr_fisherScore.RData")
# # pp <- pp2
# 
# do.call(rbind, lapply(pp, function(i) i[["clustersAddon"]][["nclusters"]])) %>% plot
# 
# # library(rlist)
# # pp_filter <- list.filter(pp, clustersAddon$nclusters < 20)
# 
# pp_filter <- pp
# 
# 
# # to get number of clusters from each Kmeans cluster
# N_clusters <- do.call(rbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["nclusters"]]))
# 
# do.call(c, lapply(pp_filter, function(i) i[["clustersAddon"]][["nclusters"]])) %>% plot
# 
# # variance explained
# do.call(c, lapply(pp_filter, function(i) i[["clustersAddon"]][["varExplained"]])) %>% plot
# 
# 
# pcTrain <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["PC"]]))
# avgTrain <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["averageExpr"]]))
# 
# pcTrain <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAll"]][["PC"]]))
# avgTrain <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAll"]][["averageExpr"]]))
# 
# pp$`1`$clustersAll$nclusters
# pp$`1`$clustersAddon$nclusters
# pp$`1`$clustersAddon$PC %>% head
# pp$`1`$clustersAll$PC %>% head
# 
# colnames(pcTrain) <- gsub("\\.","_",colnames(pcTrain))
# colnames(pcTrain) <- paste0("PC",colnames(pcTrain))
# datt <- cbind(pcTrain, Y = Y[trainIndex])
# colnames(datt)
# str(datt)
# 
# dim(pcTrain)
# pcTest <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["PCTest"]]))
# avgTest <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["averageExprTest"]]))
# dim(pcTest)
# 
# 
# library(ComplexHeatmap)
# require(circlize)
# pcTrain %>% dim
# max(pcTrain);min(pcTrain)
# apply(pcTrain, 2, max) %>% as.numeric() %>% max
# apply(pcTrain, 2, min) %>% as.numeric() %>% min
# 
# 
# cm <- colorRamp2(seq(min(pcTrain), max(pcTrain), length.out = 50), viridis(50))
# ht1 = Heatmap(t(pcTrain[which(pp[[1]][["etrain"]]==0),]), 
#               name = "E=0",
#               # col = viridis(100),
#               col = cm,
#               # column_title = "E = 0 : Age [4.8, 11.3]",
#               column_title = "Income_Level: 1-7",
#               clustering_method_rows = "average",
#               # clustering_method_columns = "average",
#               show_row_names = FALSE)
# ht2 = Heatmap(t(pcTrain[which(pp[[1]][["etrain"]]==1),]), 
#               name = "E=1",
#               # col = viridis(100), 
#               col = cm,
#               clustering_method_rows = "average",
#               # clustering_method_columns = "average",
#               # column_title = "E = 1 : Age [11.3, 18]",
#               column_title = "Income_Level: 8-10",
#               show_row_names = FALSE)
# ht1 + ht2
# 
# pheatmap::pheatmap(t(pcTrain[which(pp[[1]][["etrain"]]==0),]))
# 
# max_heat <- max(c(max(pcTrain[which(pp[[1]][["etrain"]]==0),]),max(pcTrain[which(pp[[1]][["etrain"]]==1),])))
# min_heat <- min(c(min(pcTrain[which(pp[[1]][["etrain"]]==0),]),min(pcTrain[which(pp[[1]][["etrain"]]==1),])))
# 
# pheatmap::pheatmap(t(pcTrain[which(pp[[1]][["etrain"]]==0),]),
#                    clustering_method = "average",
#                    color = viridis(100),
#                    breaks = seq(min_heat, max_heat, length.out = 101))
# pheatmap::pheatmap(t(pcTrain[which(pp[[1]][["etrain"]]==1),]),
#                    clustering_method = "average",
#                    color = viridis(100),
#                    breaks = seq(min_heat, max_heat, length.out = 101))
# 
# # ,
# #                    legend_breaks = round(seq(min(min_max_heat$V2), max(min_max_heat$V1), length.out = 12),1),
# #                    legend_labels = round(seq(min(min_max_heat$V2), max(min_max_heat$V1), length.out = 12),1))
# 
# cv.fit <- cv.glmnet(x = as.matrix(pcTrain), y = Y[trainIndex], alpha = 0.5, standardize = T, intercept=T)
# cv.fit <- cv.glmnet(x = as.matrix(avgTrain), y = Y[trainIndex], alpha = 1)
# 
# plot(cv.fit)
# as.matrix(coef(cv.fit, s = "lambda.1se"))[which(as.matrix(coef(cv.fit, s = "lambda.1se"))!=0),,drop=F]
# as.matrix(coef(cv.fit, s = "lambda.min"))[which(as.matrix(coef(cv.fit, s = "lambda.min"))!=0),,drop=F]
# cv.fit$cvm[which.min(cv.fit$cvm)]
# 
# 
# dim(pcTrain)
