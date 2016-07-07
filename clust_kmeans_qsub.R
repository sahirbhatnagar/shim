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

uni_results[, cv:=sd/mean, by = probe]

theta = seq(from=0, to=0.99, by=0.02)
filterChoices = data.frame(
  `mean` = uni_results$mean,
  `geneID` = 1:nrow(uni_results),
  `sd` = uni_results$sd,
  `cv` = uni_results$cv
)

rejChoices = sapply(filterChoices, function(f)
  filtered_R(alpha=0.1, filter=f, test=uni_results$pvalue, theta=theta, method="BH"))

# library("RColorBrewer")
# myColours = brewer.pal(ncol(filterChoices), "Set1")
# matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
#         xlab=expression(theta), ylab="number of rejections")
# legend("topleft", legend=colnames(filterChoices), fill=myColours)
# abline(v = theta[order(rejChoices[,"mean"], decreasing = TRUE)][1], lty=2)
# 
# theta[order(rejChoices[,"mean"], decreasing = TRUE)]
# quantile(uni_results$mean)

thetaThreshold <- theta[order(rejChoices[,"mean"], decreasing = TRUE)][1]
filtered_probes <- uni_results[with(uni_results, mean > quantile(mean, thetaThreshold)),]$probe
head(filtered_probes)

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

load("kclustData_NIH_50_kmeans_after_mean_filter.RData")
table(kclustData$result)

# these contain the brain probes in their respective kmeans clusters
# subjects are rows, columns are probes
clusters <- lapply(by(t(thicknessMat),kclustData$result,identity), function(i) t(as.matrix(i)))
lapply(clusters, dim)
lapply(clusters, is.matrix)

# lapply(clusters, function(i) all(colnames(i %in% brain_probes))) %>% unlist() %>% all

# Create training and test set indices ------------------------------------

Y <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, phenotypeVariable, with = F] %>% as.matrix() %>% drop

trainIndex <- caret::createDataPartition(Y, p = 1, list = FALSE, times = 1) %>% drop
testIndex <- trainIndex
# testIndex <- which(seq_len(length(Y)) %ni% trainIndex)
# all(c(trainIndex, testIndex) %in% seq_len(length(Y)))

# we also need the exposure variable as a vector
E <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, exposureVariable, with = F] %>% as.matrix() %>% drop


doMC::registerDoMC(cores = 5)
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
               mc.cores = 5)

save(pp, file = "PC_NIHPD_based_on_all_data_converted_to_50_kmeans_clusters_filtered_mean_TOM_DIFFTOM.RData")


doMC::registerDoMC(cores = 5)
pp2 <- mclapply(clusters, function(i) cluster_kmeans(data = i, 
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
                                                    nPC = 1),
               mc.cores = 5)

save(pp2, file = "PC_NIHPD_based_on_all_data_converted_to_50_kmeans_clusters_filtered_mean_corr_fisherScore.RData")


