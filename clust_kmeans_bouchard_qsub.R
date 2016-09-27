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

rm(list=ls())
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/packages.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/data_cleaning.R")

Rcpp::sourceCpp("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/kmeans.cpp")
options(scipen = 999, digits = 4)

# load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/kclustData.RData")

# these contain the cpg probes in their respective kmeans clusters
# subjects are rows, columns are probes
# clusters <- lapply(by(placentaALL,kclustData$result,identity), function(i) t(as.matrix(i)))
# lapply(clusters, dim)

probe_sd <- rowSds(placentaALL)

# 10,000 most variable probes
# filterd_probes <- probe_sd[probe_sd > quantile(probe_sd,0.956485)] %>% names
# filterd_probes %>% length()

# 5,000 most variable probes
filterd_probes <- probe_sd[probe_sd > quantile(probe_sd,0.978245)] %>% names
filterd_probes %>% length()


# parameterIndex <- as.numeric(commandArgs(trailingOnly = T))
# parameterIndex

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

# doMC::registerDoMC(cores = 20)
# pp <- mclapply(clusters[parameterIndex], function(i) cluster_kmeans(data = i, 
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
#                mc.cores = 20)


pp <- cluster_kmeans(data = t(placentaALL[filterd_probes,]), 
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


save(pp, file = paste0(
  "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/","1_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM_Sept19.RData"
))

# save(pp, file = paste0(
#   "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git",tempfile(tmpdir = ""),"_2_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM.RData"
#   )
# )
# 
# pp2 <- cluster_kmeans(data = t(placentaALL[filterd_probes,]), 
#                      exposure = E,
#                      response = Y,
#                      min_cluster_size = 50,
#                      train_index = trainIndex, 
#                      test_index = testIndex,
#                      cluster_distance = "corr",
#                      cluster_method = "hclust",
#                      cut_method = "dynamic",
#                      agglomeration_method = "average",
#                      distance_method = "euclidean",
#                      eclust_add = TRUE,
#                      eclust_distance = "fisherScore",
#                      nPC = 2)
# 
# 
# save(pp2, file = paste0(
#   "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git",tempfile(tmpdir = ""),"_2_PC_bouchard_sd_filter_10K_probes_corr_fisherScore.RData"
# )
# )



# Mean filter - 40k probes ------------------------------------------------

# save(pp, file = paste0(
#   "/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git",tempfile(tmpdir = ""),"PC_bouchard_based_on_all_data_converted_to_50_kmeans_clusters.RData"
#   )
# )
# 
# 
# rm(pl)
# load("file75862e2e3908PC_bouchard_based_on_all_data_converted_to_50_kmeans_clusters.RData")
# c1 <- pp
# rm(pp)
# load("file336f78f6664cPC_bouchard_based_on_all_data_converted_to_50_kmeans_clusters.RData")
# c2 <- pp
# rm(pp)
# pp <- c(c1,c2)
# 
# 
# pp[[1]] %>% str
# 
# 
# library(rlist)
# pp_filter <- list.filter(pp, clustersAddon$nclusters < 500)
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
# # to combine all the principal components
# pcTrain <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["PC"]]))
# avgTrain <- do.call(cbind, lapply(pp_filter, function(i) i[["clustersAddon"]][["averageExpr"]]))
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
# library(ComplexHeatmap)
# require(circlize)
# 
# cm <- colorRamp2(seq(min(pcTrain), max(pcTrain), length.out = 100), viridis(100))
# ht1 = Heatmap(t(pcTrain[which(pp_filter[[1]][["etrain"]]==0),]), 
#               name = "E=0",
#               # col = viridis(10), 
#               col = cm,
#               # column_title = "E = 0 : Age [4.8, 11.3]",
#               # column_title = "Income_Level: 1-7",
#               column_title = "NGD",
#               show_row_names = FALSE)
# ht2 = Heatmap(t(pcTrain[which(pp_filter[[1]][["etrain"]]==1),]), 
#               name = "E=1",
#               # col = viridis(10), 
#               col = cm,
#               # column_title = "E = 1 : Age [11.3, 18]",
#               # column_title = "Income_Level: 8-10",
#               column_title = "GD",
#               show_row_names = FALSE)
# ht1 + ht2
# 
# cm <- colorRamp2(seq(min(avgTrain), max(avgTrain), length.out = 100), viridis(100))
# ht1 = Heatmap(t(avgTrain[which(pp_filter[[1]][["etrain"]]==0),]), 
#               name = "E=0",
#               # col = viridis(10), 
#               col = cm,
#               # column_title = "E = 0 : Age [4.8, 11.3]",
#               # column_title = "Income_Level: 1-7",
#               column_title = "NGD",
#               show_row_names = FALSE)
# ht2 = Heatmap(t(avgTrain[which(pp_filter[[1]][["etrain"]]==1),]), 
#               name = "E=1",
#               # col = viridis(10), 
#               col = cm,
#               # column_title = "E = 1 : Age [11.3, 18]",
#               # column_title = "Income_Level: 8-10",
#               column_title = "GD",
#               show_row_names = FALSE)
# ht1 + ht2
# 
# 
# 
# 
# 
# 
# 
# 
# cv.fit <- cv.glmnet(x = as.matrix(pcTrain), y = Y[trainIndex], alpha = 0.5, standardize = T, intercept=T)
# # cv.fit <- cv.glmnet(x = as.matrix(avgTrain), y = Y[trainIndex], alpha = 0.5)
# 
# plot(cv.fit)
# as.matrix(coef(cv.fit, s = "lambda.1se"))[which(as.matrix(coef(cv.fit, s = "lambda.1se"))!=0),,drop=F]
# as.matrix(coef(cv.fit, s = "lambda.min"))[which(as.matrix(coef(cv.fit, s = "lambda.min"))!=0),,drop=F]
# cv.fit$cvm[which.min(cv.fit$cvm)]
# 
# 
# # Create Interaction Data for cv.shim ----------------------------------------
# 

# 
# kl <- prepare_data(data = cbind(pcTrain, pheno = Y[trainIndex], income = E[trainIndex]),
#                    response = "pheno", exposure = "income")
# 
# kl$X %>% dim()
# kl$main_effect_names
# kl$interaction_names
# 
# # cv.fit2 <- cv.glmnet(x = kl$X, y = kl$Y, alpha = 0.5)
# # plot(cv.fit2)
# # as.matrix(coef(cv.fit2, s = "lambda.1se"))[which(as.matrix(coef(cv.fit2, s = "lambda.1se"))!=0),,drop=F]
# # as.matrix(coef(cv.fit2, s = "lambda.min"))[which(as.matrix(coef(cv.fit2, s = "lambda.min"))!=0),,drop=F]
# 
# library(eclust)
# library(doMC)
# registerDoMC(cores = 10)
# system.time(shim <- shim(x =kl$X, y = kl$Y,
#                                main.effect.names = kl$main_effect_names,
#                                interaction.names = kl$interaction_names,
#                                verbose = TRUE))
# 
# system.time(cv_shim <- cv.shim(x =kl$X, y = kl$Y,
#                                main.effect.names = kl$main_effect_names,
#                                interaction.names = kl$interaction_names,
#                                parallel = TRUE, verbose = TRUE,
#                                type.measure = c("mse"), 
#                                nfolds = 10))
# 





# # Sd filter - 10K probes --------------------------------------------------
# 
# # plot var explained ----
# 
# # load("file1813f434eddPC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM.RData")
# # load("file18134f428062PC_bouchard_sd_filter_10K_probes_corr_fisherScore.RData")
# 
# rm(pp, pp2)
# 
# load("file97672a4c93c72_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM.RData")
# load("file97674b79bc052_PC_bouchard_sd_filter_10K_probes_corr_fisherScore.RData")
# 
# # to combine all the principal components
# pcTrain_TOM <- pp$clustersAddon$PC
# dim(pcTrain_TOM)
# head(pcTrain_TOM)
# avgTrain_TOM <- pp$clustersAddon$averageExpr
# 
# # to combine all the principal components
# pcTrain_corr <- pp2$clustersAddon$PC
# dim(pcTrain_corr)
# head(pcTrain_corr)
# avgTrain_corr <- pp2$clustersAddon$averageExpr
# 
# varexp_PC1_TOM <- pp$clustersAddon$varExplained[seq(1, length(pp$clustersAddon$varExplained), by = 2)]
# varexp_PC2_TOM <- pp$clustersAddon$varExplained[seq(2, length(pp$clustersAddon$varExplained), by = 2)]
# 
# varexp_PC1_corr <- pp2$clustersAddon$varExplained[seq(1, length(pp2$clustersAddon$varExplained), by = 2)]
# varexp_PC2_corr <- pp2$clustersAddon$varExplained[seq(2, length(pp2$clustersAddon$varExplained), by = 2)]
# 
# 
# dTOM <- data.frame(index = seq_len(length(varexp_PC1_TOM)), varexp_PC1_TOM, varexp_PC2_TOM) %>% 
#   gather(type, value, -index) %>% 
#   separate(type, c("measure", "PC", "type"))
# dcorr <- data.frame(index = seq_len(length(varexp_PC1_corr)), varexp_PC1_corr, varexp_PC2_corr) %>% 
#   gather(type, value, -index) %>% 
#   separate(type, c("measure", "PC", "type"))
# 
# var_expl_data <- rbind(dTOM, dcorr)
# 
# p <- ggplot(var_expl_data, aes(x = index, y = value, color = PC))
# p + geom_point(size=2) + facet_wrap(~type) + ylab("variance explained") + theme_bw()
# 
# # plot(seq_len(length(var_exp_pc1)), var_exp_pc1, pch = 19, col = "red", 
# #      ylim = c(0, max(var_exp_pc1, var_exp_pc2)),
# #      xlab = "cluster index", ylab = "variance explained")
# # points(seq_len(length(var_exp_pc1)), var_exp_pc2, pch = 19, col = "blue")
# # legend("topright", c("PC1", "PC2"), col = c("red","blue"),
# #        pch = c(19,19), bg = "gray90")
# 
# colnames(pcTrain_TOM) <- gsub("\\.","_",colnames(pcTrain_TOM))
# # colnames(pcTrain_TOM) <- paste0("PC",colnames(pcTrain_TOM))
# datt <- cbind(pcTrain_TOM, Y = Y[trainIndex], age = DT.pheno.placenta$Age_gestationnel, sex = DT.pheno.placenta$Sexe)
# colnames(datt)
# str(datt)
# 
# 
# 
# # sent to Celia July 13
# bouchard_2PC_10K_probes_TOM_DIFFTOM <- cbind(pcTrain_TOM, 
#                                             bmizscore = Y[trainIndex], 
#                                             gdstatus = E[trainIndex])
# 
# write.csv(bouchard_2PC_10K_probes_TOM_DIFFTOM, file = "bouchard_2PC_10K_TOM_DIFFTOM.csv", 
#           quote = FALSE,row.names = FALSE)
# head(bouchard_2PC_10K_probes_TOM_DIFFTOM)
# 
# 
# max_heat <- max(c(max(pcTrain_TOM[which(pp$etrain==0),]),max(pcTrain_TOM[which(pp$etrain==1),])))
# min_heat <- min(c(min(pcTrain_TOM[which(pp$etrain==0),]),min(pcTrain_TOM[which(pp$etrain==1),])))
# 
# pheatmap::pheatmap(t(pcTrain_TOM[which(pp$etrain==0),]),
#                    clustering_method = "average",
#                    color = viridis(100),
#                    breaks = seq(min_heat, max_heat, length.out = 101))
# pheatmap::pheatmap(t(pcTrain_TOM[which(pp$etrain==1),]),
#                    clustering_method = "average",
#                    color = viridis(100),
#                    breaks = seq(min_heat, max_heat, length.out = 101))
# 
# 
# library(ComplexHeatmap)
# require(circlize)
# 
# cm <- colorRamp2(seq(min_heat, max_heat, length.out = 100), viridis(100))
# ht1 = Heatmap(t(pcTrain_TOM[which(pp$etrain==0),]), 
#               name = "E=0",
#               # col = viridis(10), 
#               col = cm,
#               # column_title = "E = 0 : Age [4.8, 11.3]",
#               # column_title = "Income_Level: 1-7",
#               column_title = "NGD",
#               show_row_names = FALSE)
# ht2 = Heatmap(t(pcTrain_TOM[which(pp$etrain==1),]), 
#               name = "E=1",
#               # col = viridis(10), 
#               col = cm,
#               # column_title = "E = 1 : Age [11.3, 18]",
#               # column_title = "Income_Level: 8-10",
#               column_title = "GD",
#               show_row_names = FALSE)
# ht1 + ht2
# 
# 
# cv.fit <- cv.glmnet(x = as.matrix(pcTrain_TOM), y = Y[trainIndex], alpha = 0.5, standardize = T, intercept=T)
# # cv.fit <- cv.glmnet(x = as.matrix(avgTrain), y = Y[trainIndex], alpha = 0.5)
# 
# plot(cv.fit)
# as.matrix(coef(cv.fit, s = "lambda.1se"))[which(as.matrix(coef(cv.fit, s = "lambda.1se"))!=0),,drop=F]
# as.matrix(coef(cv.fit, s = "lambda.min"))[which(as.matrix(coef(cv.fit, s = "lambda.min"))!=0),,drop=F]
# 
# kl <- prepare_data(data = cbind(pcTrain_TOM, pheno = Y[trainIndex], 
#                                 income = E[trainIndex]
# #                                 ,
# #                                 age = DT.pheno.placenta$Age_gestationnel, 
# #                                 sex = DT.pheno.placenta$Sexe
#                                 ),
#                    response = "pheno", exposure = "income")
# 
# kl$X %>% dim()
# kl$main_effect_names
# kl$interaction_names
# 
# library(eclust)
# library(doMC)
# registerDoMC(cores = 10)
# system.time(shim <- shim(x =kl$X, y = kl$Y,
#                                main.effect.names = kl$main_effect_names,
#                                interaction.names = kl$interaction_names,
#                                verbose = TRUE))
# 
# system.time(cv_shim <- cv.shim(x =kl$X, y = kl$Y,
#                                main.effect.names = kl$main_effect_names,
#                                interaction.names = kl$interaction_names,
#                                parallel = TRUE, verbose = TRUE,
#                                type.measure = c("mse"), 
#                                nfolds = 10))
# plot(cv_shim)
# coef(cv_shim)
# cv_shim$shim.fit
# 
# 
# cv.fit2 <- cv.glmnet(x = kl$X, y = kl$Y, alpha = 0.5)
# plot(cv.fit2)
# as.matrix(coef(cv.fit2, s = "lambda.1se"))[which(as.matrix(coef(cv.fit2, s = "lambda.1se"))!=0),,drop=F]
# as.matrix(coef(cv.fit2, s = "lambda.min"))[which(as.matrix(coef(cv.fit2, s = "lambda.min"))!=0),,drop=F]
# 


# MARS --------------------------------------------------------------------

# this loads an object called pp
# load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/file19e17f13c434_2_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM.RData")


load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/1_PC_bouchard_sd_filter_5K_probes_TOM_DIFFTOM_Sept19.RData")
res5k <- pp

rm(pp)

load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/1_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM_Sept19.RData")

res10k <- pp
rm(pp)


load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/1_PC_bouchard_sd_filter_20K_probes_TOM_DIFFTOM_Sept19.RData")
res20k <- pp
rm(pp)


# this loads an object called pp2
# load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/file19e150ebfb6e_2_PC_bouchard_sd_filter_10K_probes_corr_fisherScore.RData")

# to combine all the principal components

# 5k most variable
pcTrain_TOM <- res5k$clustersAddon$PC
res5k$clustersAll$nclusters
res5k$clustersAddon$nclusters

res5k$clustersAddon$prcompObj %>% str

dim(pcTrain_TOM)
head(pcTrain_TOM)
avgTrain_TOM <- res5k$clustersAddon$averageExpr


pcTrain_TOM_clust <- res5k$clustersAll$PC
res5k$clustersAll$nclusters
res5k$clustersAddon$nclusters

res5k$clustersAll$prcompObj %>% str

dim(pcTrain_TOM_clust)
head(pcTrain_TOM_clust)
avgTrain_TOM_clust <- res5k$clustersAll$averageExpr




# 10k most variable
pcTrain_TOM_10k <- res10k$clustersAddon$PC
res10k$clustersAll$nclusters
res10k$clustersAddon$nclusters

res10k$clustersAddon$prcompObj %>% str

dim(pcTrain_TOM_10k)
head(pcTrain_TOM_10k)
avgTrain_TOM_10k <- res10k$clustersAddon$averageExpr

# 20k most variable
pcTrain_TOM_20k <- res20k$clustersAddon$PC
res20k$clustersAll$nclusters
res20k$clustersAddon$nclusters

res20k$clustersAddon$prcompObj %>% str

dim(pcTrain_TOM_20k)
head(pcTrain_TOM_20k)
avgTrain_TOM_20k <- res20k$clustersAddon$averageExpr


# to combine all the principal components
# pcTrain_corr <- pp2$clustersAddon$PC
# dim(pcTrain_corr)
# head(pcTrain_corr)
# avgTrain_corr <- pp2$clustersAddon$averageExpr
# 
# varexp_PC1_TOM <- pp$clustersAddon$varExplained[seq(1, length(pp$clustersAddon$varExplained), by = 2)]
# varexp_PC2_TOM <- pp$clustersAddon$varExplained[seq(2, length(pp$clustersAddon$varExplained), by = 2)]
# 
# varexp_PC1_corr <- pp2$clustersAddon$varExplained[seq(1, length(pp2$clustersAddon$varExplained), by = 2)]
# varexp_PC2_corr <- pp2$clustersAddon$varExplained[seq(2, length(pp2$clustersAddon$varExplained), by = 2)]
# 
# 
# dTOM <- data.frame(index = seq_len(length(varexp_PC1_TOM)), varexp_PC1_TOM, varexp_PC2_TOM) %>% 
#   gather(type, value, -index) %>% 
#   separate(type, c("measure", "PC", "type"))
# dcorr <- data.frame(index = seq_len(length(varexp_PC1_corr)), varexp_PC1_corr, varexp_PC2_corr) %>% 
#   gather(type, value, -index) %>% 
#   separate(type, c("measure", "PC", "type"))
# 
# var_expl_data <- rbind(dTOM, dcorr)
# 
# p <- ggplot(var_expl_data, aes(x = index, y = value, color = PC))
# p + geom_point(size=2) + facet_wrap(~type) + ylab("variance explained") + theme_bw()

max_heat <- max(c(max(pcTrain_TOM[which(res5k$etrain==0),]),max(pcTrain_TOM[which(res5k$etrain==1),])))
min_heat <- min(c(min(pcTrain_TOM[which(res5k$etrain==0),]),min(pcTrain_TOM[which(res5k$etrain==1),])))

library(ComplexHeatmap)
require(circlize)

cm <- colorRamp2(seq(min_heat, max_heat, length.out = 100), viridis(100))
ht1 = Heatmap(t(pcTrain_TOM[which(res5k$etrain==0),]), 
              name = "E=0",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 0 : Age [4.8, 11.3]",
              # column_title = "Income_Level: 1-7",
              column_title = "NGD",
              show_row_names = FALSE)
ht2 = Heatmap(t(pcTrain_TOM[which(res5k$etrain==1),]), 
              name = "E=1",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 1 : Age [11.3, 18]",
              # column_title = "Income_Level: 8-10",
              column_title = "GD",
              show_row_names = FALSE)

png(paste0("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/","1_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM_Sept19.png"))
ht1 + ht2
dev.off()



max_heat <- max(c(max(pcTrain_TOM_10k[which(res10k$etrain==0),]),max(pcTrain_TOM_10k[which(res10k$etrain==1),])))
min_heat <- min(c(min(pcTrain_TOM_10k[which(res10k$etrain==0),]),min(pcTrain_TOM_10k[which(res10k$etrain==1),])))
cm <- colorRamp2(seq(min_heat, max_heat, length.out = 100), viridis(100))
ht1 = Heatmap(t(pcTrain_TOM_10k[which(res10k$etrain==0),]), 
              name = "E=0",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 0 : Age [4.8, 11.3]",
              # column_title = "Income_Level: 1-7",
              column_title = "NGD",
              show_row_names = FALSE)
ht2 = Heatmap(t(pcTrain_TOM_10k[which(res10k$etrain==1),]), 
              name = "E=1",
              # col = viridis(10), 
              col = cm,
              # column_title = "E = 1 : Age [11.3, 18]",
              # column_title = "Income_Level: 8-10",
              column_title = "GD",
              show_row_names = FALSE)

# png(paste0("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/","1_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM_Sept19.png"))
ht1 + ht2
# dev.off()



kl_5k <- prepare_data(data = cbind(pcTrain_TOM, pheno = Y[trainIndex],
                                income = E[trainIndex]),
                   response = "pheno", exposure = "income")

kl_10k <- prepare_data(data = cbind(pcTrain_TOM_10k, pheno = Y[trainIndex],
                                   income = E[trainIndex]),
                      response = "pheno", exposure = "income")

kl_20k <- prepare_data(data = cbind(pcTrain_TOM_20k, pheno = Y[trainIndex],
                                    income = E[trainIndex]),
                       response = "pheno", exposure = "income")


kl_5k$X %>% dim
kl_5k$X[,1:45] %>% head

kl_10k$X %>% dim
kl_10k$X[,1:78] %>% head

kl_20k$X %>% dim
kl_20k$X[,1:142] %>% head


fit.earth.5k <- earth::earth(x = kl_5k$X[,1:45], y = kl_5k$Y, pmethod = "backward", keepxy = TRUE, 
                             degree = 3, nfold = 5, ncross = 3, trace = 4, nk = 1000)
summary(fit.earth.5k)
plot(fit.earth.5k)
plotmo(fit.earth.5k)
evimp(fit.earth.5k)





# Performance using Cross Validation --------------------------------------

fitControl <-  trainControl(method = "repeatedcv",
                            number = 10,
                            repeats = 5,
                            verboseIter = TRUE)

# fitControl <-  trainControl(method = "cv",
#                             number = 5,
#                             repeats = 1,
#                             verboseIter = TRUE)


# Define the candidate models to test
marsGrid <- expand.grid(.degree = 1:3, .nprune = 1000)
# Fix the seed so that the results can be reproduced

# MARS tuned ECLUST ----
set.seed(1056)
marsEclust <- train(kl_5k$X[,1:45], kl_5k$Y,
                   method = "earth",
                   trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                   # Explicitly declare the candidate models to test
                   tuneGrid = marsGrid,
                   trControl = fitControl)
# marsTuned$finalModel$bx
# plot(marsTuned)
# marsTuned %>% summary()
# marsTuned$bestTune
# 
# fit.earth <- earth::earth(x = kl_5k$X[,1:45], y = kl_5k$Y,
#                           trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward", 
#                           degree = marsTuned$bestTune$degree)
# summary(fit.earth)
# plotmo(fit.earth)

# MARS tuned Original ----
set.seed(1056)
marsOriginal <- train(t(placentaALL[filterd_probes,]), kl_5k$Y,
                   method = "earth",
                   trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                   # Explicitly declare the candidate models to test
                   tuneGrid = marsGrid,
                   trControl = fitControl)

# MARS tuned bagged ECLUST ----
set.seed(1056)
bagMarsEclust <- train(kl_5k$X[,1:45], kl_5k$Y,
                   method = "bagEarth", B = 10,
                   trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                   # Explicitly declare the candidate models to test
                   tuneGrid = marsGrid,
                   trControl = fitControl)

# marsTunedBag$finalModel$fit$Resample10
# summary(marsTunedBag)
# 
# predict(marsTunedBag)
# predict(marsTuned)


# Lasso + ECLUST ----
# apparently this does tune over both alpha and lambda, but just need to provide lambda
# http://stats.stackexchange.com/questions/69638/does-caret-train-function-for-glmnet-cross-validate-for-both-alpha-and-lambda?rq=1
# lassoGrid <- expand.grid(.alpha = seq(0,1,length.out = 5), .lambda = seq(min(glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda),
#                                                               max(glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda),length.out = 5))
# glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda

set.seed(1056)
lassoEclust <- train(x = kl_5k$X, y = kl_5k$Y,
                  method = "glmnet",
                  tuneLength = 20,
                  trControl = fitControl)

# Lasso Original ----

set.seed(1056)
lassoOriginal <- train(t(placentaALL[filterd_probes,]), kl_5k$Y,
                     method = "glmnet",
                     tuneLength = 20,
                     trControl = fitControl)





resamps <- resamples(list(earth = marsTuned,
                          bagEarth = marsTunedBag))


trellis.par.set("theme1")
bwplot(resamps, layout = c(3, 1))

summary(resamps)

resamps


data(trees)
fit1 <- bagEarth(Volume ~ ., data = trees, keepX = TRUE)
fit2 <- bagEarth(Volume ~ ., data = trees, keepX = FALSE)
hist(predict(fit1) - predict(fit2))







# MARS 5k Original ----

t(placentaALL[filterd_probes,]) %>% str
fit.earth.5k_original <- earth::earth(x = t(placentaALL[filterd_probes,]), y = Y, pmethod = "backward", keepxy = TRUE, 
                                      degree = 3, trace = 4, nk = 1000, nfold = 5, ncross = 3)
summary(fit.earth.5k_original)
plot(fit.earth.5k_original)
plotmo(fit.earth.5k_original)
evimp(fit.earth.5k_original)




fit.earth.10k <- earth::earth(x = kl_10k$X[,1:78], y = kl_10k$Y, pmethod = "backward", keepxy = TRUE, degree = 3, nfold = 5, ncross = 3, trace = 4)
summary(fit.earth.10k)
plot(fit.earth.10k)
plotmo(fit.earth.10k)
evimp(fit.earth.10k)


fit.earth.10k_original <- earth::earth(x = t(placentaALL[filterd_probes,]), y = Y, pmethod = "backward", keepxy = TRUE, 
                                      degree = 3, trace = 4, nk = 1000, nfold = 3, ncross = 3)
summary(fit.earth.10k_original)
plot(fit.earth.10k_original)
plotmo(fit.earth.10k_original)
evimp(fit.earth.10k_original)



fit.earth.20k <- earth::earth(x = kl_20k$X[,1:142], y = kl_20k$Y, pmethod = "backward", keepxy = TRUE, degree = 3, nfold = 5, ncross = 3, trace = 4)
summary(fit.earth.20k)
plot(fit.earth.20k)
plotmo(fit.earth.20k)
evimp(fit.earth.20k)


pc_21_5k <- rownames(res5k$clustersAddon$prcompObj[[21]]$rotation)
pc_27_5k <- rownames(res5k$clustersAddon$prcompObj[[27]]$rotation)

pc_56_10k <- rownames(res10k$clustersAddon$prcompObj[[56]]$rotation)

pc_111_20k <- rownames(res20k$clustersAddon$prcompObj[[111]]$rotation)
pc_114_20k <- rownames(res20k$clustersAddon$prcompObj[[114]]$rotation)
pc_128_20k <- rownames(res20k$clustersAddon$prcompObj[[128]]$rotation)


genes_5k <- DT.placenta[rn %in% c(pc_21_5k,pc_27_5k)]$nearestGeneSymbol %>% unique
genes_5k %>% length()

genes_10k <- DT.placenta[rn %in% c(pc_56_10k)]$nearestGeneSymbol %>% unique
genes_10k %>% length()

genes_20k <- DT.placenta[rn %in% c(pc_111_20k, pc_114_20k, pc_128_20k)]$nearestGeneSymbol %>% unique
genes_20k %>% length()


library(VennDiagram)
overlap <- VennDiagram::calculate.overlap(x = list('5k' = genes_5k, 
                                    '10k' = genes_10k,
                                    '20k' = genes_20k))


VennDiagram::venn.diagram(list('5k' = genes_5k, 
                               '10k' = genes_10k,
                               '20k' = genes_20k), filename = "overlap.png", imagetype = "png")

intersect(genes_5k,genes_10k) %>% unique()




# enter in gene mania

write.table(genes_10k,
quote = F, sep = "\n", row.names = F, col.names = FALSE)

# hello world


## try http:// if https:// URLs are not supported
# source("http://bioconductor.org/biocLite.R")
# biocLite("gage")













# resamples()
# 
# head(predict(marsTuned, kl$X))
# varImp(marsTuned)
# 
# earthPred <- predict(marsTuned, newdata = kl$X)
## The function ' postResample ' can be used to get the test set
## perforamnce values
# postResample(pred = earthPred, obs = kl$X)


# fitControl <-  trainControl(method = "repeatedcv", 
#                             number = 10, 
#                             repeats = 3,
#                             verboseIter = TRUE)
# 
# # Define the candidate models to test
# marsGrid <- expand.grid(.degree = 1:2, .nprune = 2:20)
# # Fix the seed so that the results can be reproduced
# set.seed(1056)
# marsTuned <- train(kl$X[,1:154], kl$Y,
#                    method = "earth",
#                    # Explicitly declare the candidate models to test
#                    tuneGrid = marsGrid,
#                    trControl = fitControl)
# marsTuned$finalModel$bx
# plot(marsTuned)
# plotmo(marsTuned)
# 
# marsTuned$bestTune
# fit.earth <- earth::earth(x = kl$X[,1:154], y = kl$Y, 
#                           pmethod = "forward", keepxy = TRUE, degree = 2, nfold = 0, trace = 4)
# plotmo(fit.earth)
# 
# set.seed(1056)
# earthFit <- train(x = kl$X[,1:154], y = kl$Y,
#                   method = "earth",
#                   tuneLength = 10,
#                   tuneGrid = marsGrid,
#                   trControl = fitControl)
# 
# summary(earthFit)
# plot(earthFit)
# earthFit$results
# 
# set.seed(1056)
# bagearthFit <- train(x = kl$X[,1:154], y = kl$Y,
#                   method = "bagEarth",
#                   tuneLength = 10,
#                   trControl = fitControl)
# 
# plot(bagearthFit)
# 
# set.seed(1056)
# lassoFit <- train(x = kl$X[,1:154], y = kl$Y,
#                   method = "glmnet",
#                   tuneLength = 10,
#                   trControl = fitControl)
# 
# resamps <- resamples(list(MARS = earthFit,
#                           bagMARS = bagearthFit,
#                           GLMNET = lassoFit))
# 
# set.seed(1056)
# relaxoFit <- train(x = kl$X[,1:154], y = kl$Y,
#                   method = "relaxo",
#                   tuneLength = 10,
#                   trControl = fitControl)
# 
# set.seed(1056)
# rbfSVM <- train(x = kl$X[,1:154], y = kl$Y,
#                   method = "svmRadial",
#                   tuneLength = 10,
#                   trControl = fitControl)
# 
# lassoFit
# resamps <- resamples(list(MARS = earthFit,
#                           bagMARS = bagearthFit,
#                           GLMNET = lassoFit,
#                           Relaxo = relaxoFit,
#                           # Ridge = ridgeFit,
#                           SVMrbf = rbfSVM))
# 
# bwplot(resamps)
# summary(resamps)
# dotplot(resamps)
# splom(resamps)
# modelDifferences <- diff(resamps)
# ?xyplot.resamples
# 
# 
# predictors(bagearthFit)
# predictors(earthFit)
# predictors(lassoFit)
# predictors(rbfSVM)
# 
# trellis.par.set(caretTheme())
# plot(bagearthFit, type = c("g", "o"))
# 
# 
# plot(earthFit)
# summary(earthFit)
# earthFit$finalModel
# varImp(earthFit)
# str(earthFit)
# 
# # or load pre-calculated results using:
# load(url("http://caret.r-forge.r-project.org/exampleModels.RData"))
# 
# resamps <- resamples(list(CART = rpartFit,
#                           CondInfTree = ctreeFit,
#                           MARS = earthFit))
# 
# bwplot(resamps)
# summary(resamps)
# dotplot(resamps)
