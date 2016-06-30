##################################
# R source code file for analyzing NIHPD data
# Created by Sahir, May 24
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

# brain_probes[1:10]
# DT_with_pheno[,"Subject_ID", with = F]
# colnames(DT_with_pheno)[1:40]
# DT_with_pheno[,"Subject_Gender",with=F] %>% str
# DT_with_pheno[is.na(income)] %>% dim
# DT_with_pheno[is.na(WASI_Full_Scale_IQ)] %>% dim
# DT_with_pheno[is.na(Subject_Gender)] %>% dim
# DT_with_pheno[is.na(Site_Location)] %>% dim
# DT_with_pheno[!is.na(E2)] %>% dim


## ---- test-stat ----
# registerDoMC(cores = 30)
# res <- mclapply(brain_probes, function(i) {
#   
#   # i = brain_probes[1]
#   dat <- DT_with_pheno[,c("WASI_Full_Scale_IQ","Subject_Gender","Site_Location",i),with=F]
#   dat <- dat[!is.na(WASI_Full_Scale_IQ)]
#   include_E <- FALSE
#   include_interaction <- FALSE
#   fit <- lm(as.formula(paste("WASI_Full_Scale_IQ ~",i, if (include_E) "+E", if (include_interaction) paste0("+",i,":E"),
#                              "+ Subject_Gender + Site_Location")),
#             data = dat)
#   #fit %>% summary
#   coef.index <- if (include_interaction) paste0(i,":E") else i
#   
#   data.frame(probe = i,
#              pvalue = as.numeric(pt(abs(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
#                                     df = fit$df.residual, lower.tail = F)*2),
#              test.stat = as.numeric(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
#              'mean' = mean(dat[,i, with = F][[1]]),
#              'sd' = sd(dat[,i,with = F][[1]]), stringsAsFactors = FALSE)
# }, mc.cores = 30)
# 
# uni_results <- do.call(rbind, res)
# uni_results <- as.data.table(uni_results)
# save(uni_results, file = "uni_results.RData")
load("uni_results.RData")



# uni_results[order(pvalue)[1:1000]]
# uni_results[, hist(sd)]
# uni_results[, hist(test.stat)]
# str(DT)
# DT <- as.data.table(dat$Thick.81924)
# DT[, lapply(.SD, cv)] %>% as.numeric() %>%  hist() 

## ---- filter-stats ----

# with(uni_results,
#      plot(rank(mean)/length(mean), -log10(pvalue), pch=16, cex=0.45))
# 
# par(mfrow = c(1,2))
# with(uni_results,
#      plot(rank(sd)/length(sd), -log10(pvalue), pch=16, cex=0.45))
# 
# trsf = function(n) log10(n+1)
# plot(ecdf(trsf(uni_results$sd)), xlab=body(trsf), main="")
# dev.off()
# 
# 
# theta = seq(from=0, to=0.5, by=0.1)
# pBH = filtered_p(filter=uni_results$mean, test=uni_results$pvalue, theta=theta, method="BH")
# head(pBH)
# 
# rejection_plot(pBH, at="sample",
#                xlim=c(0, 0.5), ylim=c(0, 15000),
#                xlab="FDR cutoff (Benjamini & Hochberg adjusted p-value)", main="")
# 
# theta = seq(from=0, to=0.8, by=0.02)
# rejBH = filtered_R(alpha = 0.1,
#                    filter=uni_results$mean, 
#                    test=uni_results$pvalue, 
#                    theta=theta, method="BH")
# 
# plot(theta, rejBH, type="l",
#      xlab=expression(theta), ylab="number of rejections")


theta = seq(from=0, to=0.99, by=0.02)
filterChoices = data.frame(
  `mean` = uni_results$mean,
  `geneID` = 1:nrow(uni_results),
  `sd` = uni_results$sd
)
rejChoices = sapply(filterChoices, function(f)
  filtered_R(alpha=0.1, filter=f, test=uni_results$pvalue, theta=theta, method="BH"))

library("RColorBrewer")
myColours = brewer.pal(ncol(filterChoices), "Set1")
matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
        xlab=expression(theta), ylab="number of rejections")
legend("topleft", legend=colnames(filterChoices), fill=myColours)


theta[order(rejChoices[,"sd"], decreasing = TRUE)]
quantile(uni_results$sd)
thetaThreshold = theta[order(rejChoices[,"sd"], decreasing = TRUE)][2]
filtered_probes <- uni_results[with(uni_results, mean > quantile(mean, 0.939)),]$probe
filtered_probes <- uni_results[with(uni_results, sd > quantile(sd, 0.939)),]$probe

# h1 = hist(uni_results$pvalue[!pass], breaks=50, plot=FALSE)
# h2 = hist(uni_results$pvalue[pass], breaks=50, plot=FALSE)
# colori <- c(`do not pass`="khaki", `pass`="powderblue")
# 
# barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
#         col = colori, space = 0, main = "", ylab="frequency")
# text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
#      adj = c(0.5,1.7), xpd=NA)
# legend("topright", fill=rev(colori), legend=rev(names(colori)))
# 
# 
# resFilt = uni_results[pass,]
# orderInPlot = order(resFilt$pvalue)
# showInPlot = (resFilt$pvalue[orderInPlot] <= 0.06)
# alpha = 0.1
# 
# plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
#      pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]),
#      xlim=c(0,4000))
# abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)

## ---- bootSVD ----

# devtools::install_github('aaronjfisher/bootSVD')
# library(bootSVD)
# library(RColorBrewer)
# 
# pc <- fastSVD(dat$Thick.81924)
# npc <- 5
# myColours = rev(brewer.pal(npc, "Reds"))
# matplot(pc$v[, 1:npc], type = 'l', lty = 1, col = myColours) #PCs from simulated data
# legend("topleft", legend = paste0("PC", 1:npc), fill = myColours)
# pc$v %>% dim
# kclustPC <- mlKmeans(t(pc$v), 40)
# WGCNA::randIndex(table(kclustPC$result, kclustData$result))

## ---- KS-filter ----


source("https://raw.githubusercontent.com/sahirbhatnagar/eclust/master/k_filter.R")

X <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ)][,brain_probes, with = F]
Y <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ)][,"WASI_Full_Scale_IQ", with = F]

colnames(DT_with_pheno)[1:30]

obj <- k.filter(x = X[,101:200, with=F],
                y = Y$WASI_Full_Scale_IQ, response.type = "continuous", method = "fused")

ks.test()
obj$k.rank
k.filter.single

warnings()

order(obj$k.rank)[1:5]

dat <- data.frame(x=t(DT.placenta.all[order(obj$k.rank)[1:5],,drop=F]), y=DT.pheno.placenta$`imc z-score`)

lm(y ~ ., dat) %>% summary()

## ---- data-prep-for-analysis ----

# phenotypeVariable: WASI_Full_Scale_IQ, WASI_Performance_IQ, WASI_Verbal_IQ
# exposureVariable: E, E2, income_binary, income_binary2
# E: age binary, 
# E2: is only the lower and upper quartile of age as defined by Budha (Cereb cortex),
# the two middle groups are NA
# income_binary: see data_cleaning.R
# income_binary2: use this one, its more balanced
phenotypeVariable <- "WASI_Full_Scale_IQ" 
exposureVariable <- "E"

# this is the cortical thickness data for the 1st timepoint only
# all the data for all timepoints is in dat$Thick.81924
# DT_with_pheno[!is.na(phenotypeVariable)][!is.na(exposureVariable)] %>% dim
thicknessMat <- DT_with_pheno[!is.na(get(phenotypeVariable))][!is.na(get(exposureVariable))][, brain_probes, with = F] %>% as.matrix()
dim(thicknessMat)
str(thicknessMat)
kclustData <- mlKmeans(thicknessMat, 50)

# c1 <- mlKmeans(thicknessMat, 50)
# c2 <- mlKmeans(thicknessMat, 50)
# WGCNA::randIndex(table(c1$result, c2$result))

table(kclustData$result)

# these contain the brain probes in their respective kmeans clusters
# subjects are rows, columns are probes
clusters <- lapply(by(t(thicknessMat),kclustData$result,identity), function(i) t(as.matrix(i)))
lapply(clusters, dim)

# clusters[[49]] %>% dimnames
# brain_probes
# plot(kclustData$result)
# DT_with_pheno[!is.na(phenotypeVariable)][, table(get(exposureVariable), useNA = "always")]

idAll <- as.character(DT_with_pheno[!is.na(phenotypeVariable)][["Subject_ID"]])
idExposed <- as.character(DT_with_pheno[!is.na(phenotypeVariable)][get(exposureVariable)==1][["Subject_ID"]])
idUnexposed <- as.character(DT_with_pheno[!is.na(phenotypeVariable)][get(exposureVariable)==0][["Subject_ID"]])

length(idExposed)
length(idUnexposed)
all(c(idExposed,idUnexposed) %in% idAll)
length(unique(c(idExposed,idUnexposed)))==length(idAll)
# filtered_probes <- uni_results[order(test.stat,decreasing = T)[1:500]]$probe

# X <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ), 
#                    c("Subject_ID","age_binary","Subject_Gender","Site_Location",brain_probes), with = F]

# this contains E and Y and the probes
DT <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ), c("WASI_Full_Scale_IQ","E",filtered_probes), with = F]
setnames(DT, "WASI_Full_Scale_IQ", "Y")
# colnames(DT)

# this contains E and the probes only
X <- DT[, c("E",filtered_probes), with = F]
dim(X)
# dimnames(X)[[1]] <- idAll
# head(X)
# str(X)

X0 <- X[E==0]
dim(X0)
# dimnames(X0)[[1]] <- idUnexposed

X1 <- X[E==1]
dim(X1)
# dimnames(X1)[[1]] <- idExposed


# TOM_all <- TOMsimilarityFromExpr(X[,filtered_probes,with=F], nThreads = 10)
# TOM_exposed <- TOMsimilarityFromExpr(X1[,filtered_probes,with=F], nThreads = 35)
# TOM_unexposed <- TOMsimilarityFromExpr(X0[,filtered_probes,with=F], nThreads = 35)
# TOM_diff <- abs(TOM_exposed-TOM_unexposed)

Y <- as.matrix(DT[,"Y",with=F])
dimnames(Y)[[1]] <- idAll
str(Y)
is.matrix(Y)

# partition-data
trainIndex <- caret::createDataPartition(Y[,"Y"], p = 2/3, 
                                         list = FALSE, times = 1)

DT_train <- as.matrix(DT[trainIndex[,1],,])
colnames(DT_train)

DT_test <- as.matrix(DT[-trainIndex[,1],,])
colnames(DT_test)

# X_train and X_test contain the environment variable
X_train <- DT_train[ , c("E",filtered_probes)]
dim(X_train)

Y_train <- DT_train[ , "Y"]
str(Y_train)
dim(Y_train)


# fit.cv <- cv.glmnet(as.matrix(X_train), Y_train)
# plot(fit.cv)

X_test <- DT_test[ , c("E",filtered_probes)]
dim(X_test)
Y_test <- DT_test[ , "Y"]

mse_null <- crossprod(mean(Y_test) - Y_test)/length(Y_test)

# training gene expression data, this does not include the environment variable
genes_e0 <- X_train[which(X_train[,"E"] == 0), filtered_probes]
dim(genes_e0)
genes_e1 <- X_train[which(X_train[,"E"] == 1), filtered_probes]
dim(genes_e1)
genes_all <- rbind(genes_e0,genes_e1)
dim(genes_all)

message("Calculating similarity matrices")

# test set gene expression data
genes_all_test <- X_test[,filtered_probes]
dim(genes_all_test)

# corr_train_e0 <- WGCNA::cor(genes_e0)
# corr_train_e1 <- WGCNA::cor(genes_e1)
# corr_train_diff <- abs(corr_train_e1 - corr_train_e0)
corr_train_all <- WGCNA::cor(genes_all)

tom_train_e0 <- WGCNA::TOMsimilarityFromExpr(genes_e0)
dim(tom_train_e0)
dimnames(tom_train_e0)[[1]] <- dimnames(corr_train_all)[[1]]
dimnames(tom_train_e0)[[2]] <- dimnames(corr_train_all)[[2]]
head(tom_train_e0)

tom_train_e1 <- WGCNA::TOMsimilarityFromExpr(genes_e1)
dimnames(tom_train_e1)[[1]] <- dimnames(corr_train_all)[[1]]
dimnames(tom_train_e1)[[2]] <- dimnames(corr_train_all)[[2]]
head(tom_train_e1)

tom_train_diff <- abs(tom_train_e1 - tom_train_e0)
dimnames(tom_train_diff)[[1]] <- dimnames(corr_train_all)[[1]]
dimnames(tom_train_diff)[[2]] <- dimnames(corr_train_all)[[2]]
head(tom_train_diff)

tom_train_all <- WGCNA::TOMsimilarityFromExpr(genes_all)
dimnames(tom_train_all)[[1]] <- dimnames(corr_train_all)[[1]]
dimnames(tom_train_all)[[2]] <- dimnames(corr_train_all)[[2]]
head(tom_train_all)
dim(tom_train_all)

message("Creating CV folds from training data")

# Folds for Cross validation
folds_train <- caret::createFolds(Y_train, k = 10, list = T)
DT_train_folds <- lapply(folds_train, function(i) DT_train[-i,])
# these have the environment variable
X_train_folds <- lapply(DT_train_folds, function(i) i[,-grep("Y",colnames(i))])
Y_train_folds <- lapply(DT_train_folds, function(i) i[,"Y"])


cluster_distance <- "tom"
clustMethod <- "hclust"
cutMethod <- "dynamic"
method <- "average"
distanceMethod <- "euclidean"
EclustAdd <- TRUE 
EclustAddDistance <- "difftom"
cutMethod <- "dynamic" #,"gap", "fixed")

message(sprintf("Calculating number of clusters based on %s using %s with %s
                linkage and the %s to determine the number of clusters",
                cluster_distance, clustMethod, method, cutMethod))

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
# the only difference here is the distanceMethod arg
res <- if (cluster_distance %in% c("diffcorr","difftom",
                                   "corScor", "tomScor","fisherScore")) {
  clusterSimilarity(x = similarity,
                    expr = genes_all,
                    exprTest = genes_all_test,
                    distanceMethod = distanceMethod,
                    clustMethod = clustMethod,
                    cutMethod = cutMethod,
                    method = method,
                    K.max = K.max, B = B, nClusters = nClusters)
} else {
  clusterSimilarity(x = similarity,
                    expr = genes_all,
                    exprTest = genes_all_test,
                    clustMethod = clustMethod,
                    cutMethod = cutMethod,
                    method = method,
                    K.max = K.max, B = B, nClusters = nClusters)
}

message(paste("Calculating number of environment clusters based on ",
              EclustAddDistance))

# clusters based on EclustAddDistance
similarityEclust <- switch(EclustAddDistance,
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


resEclust <- if (EclustAddDistance %in% c("diffcorr","difftom",
                                          "corScor", "tomScor","fisherScore")) {
  clusterSimilarity(x = similarityEclust,
                    expr = genes_all,
                    exprTest = genes_all_test,
                    distanceMethod = distanceMethod,
                    clustMethod = clustMethod,
                    cutMethod = cutMethod,
                    method = method,
                    K.max = K.max, B = B, nClusters = nClusters)
} else {
  clusterSimilarity(x = similarityEclust,
                    expr = genes_all,
                    exprTest = genes_all_test,
                    clustMethod = clustMethod,
                    cutMethod = cutMethod,
                    method = method,
                    K.max = K.max, B = B, nClusters = nClusters)
}


# we need to combine the cluster information here
# this is based on cluster_distance only
clustersAll <- copy(res$clusters)
n_clusters_All <- res$pcInfo$nclusters

message(sprintf("There are %d clusters derived from the %s similarity matrix",
                n_clusters_All, cluster_distance))

# this is based on EclustAddDistance only
n_clusters_Eclust <- resEclust$pcInfo$nclusters
clustersEclust <- copy(resEclust$clusters)

message(sprintf("There are %d clusters derived from the %s environment similarity matrix",
                n_clusters_Eclust, EclustAddDistance))

# this is based on both
n_clusters_Addon <- n_clusters_All + n_clusters_Eclust

message(sprintf("There are a total of %d clusters derived from the %s
                similarity matrix and the %s environment similarity matrix",
                n_clusters_Addon,cluster_distance,EclustAddDistance))

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
# and the clusters from the EclustAddDistance (e.g. fisherScore)
clustersAddon <- rbindlist(list(clustersAll, clustersEclust))
clustersAddon[, table(cluster)]


result <- list(similarity = similarity,
               similarityEclust = similarityEclust,
               Y = Y, 
               X0 = X0, 
               X1 = X1, 
               X_train = X_train, 
               X_test = X_test,
               Y_train = Y_train, 
               Y_test = Y_test, 
               DT_train = DT_train,
               DT_test = DT_test, 
               n_clusters_All = n_clusters_All,
               n_clusters_Eclust = n_clusters_Eclust,
               n_clusters_Addon = n_clusters_Addon,
               clustersAll = clustersAll,
               clustersAddon = clustersAddon,
               clustersEclust = clustersEclust,
               tom_train_all = tom_train_all, 
               tom_train_diff = tom_train_diff,
               tom_train_e1 = tom_train_e1,
               tom_train_e0 = tom_train_e0,
               corr_train_all = corr_train_all,
               mse_null = mse_null, 
               DT_train_folds = DT_train_folds,
               X_train_folds = X_train_folds, 
               Y_train_folds = Y_train_folds,
               S0 = colnames(X_train)[1:5],
               beta_truth = c("dhskhd"=1.5))




## ---- univariate-p-value ----

# not working, dont use, argue that this method is ridiculous to even compare
# because there are too many choices to be made
# p <- ncol(result[["X_train"]])
# percent <- if (p > 1000) 0.005 else 0.01
# uni_res <- uniFitRda(train = result[["DT_train"]], 
#                      test = result[["DT_test"]],
#                      covariate_names = filtered_probes, # this should not include E
#                      percent = percent, stability = F,
#                      include_E = T,
#                      include_interaction = T,
#                      filter_var = F, 
#                      p = p,
#                      s0 = result[["S0"]], 
#                      true_beta = result[["beta_truth"]])


## ---- cluster-and-regress----

# we will treat the clusters as fixed i.e., even if we filter, or
# do cross validation, the group labels are predetermined by the
# above clustering procedure
# This method is based on clusters derived without accounting for the
# environment

print("starting cluster and regress with interaction")
p <- ncol(result[["X_train"]]) -  1 
clust_res <- mapply(clust_fun,
                    #summary = rep(c("pc","spc","avg"), each = 3),
                    #model = rep(c("lm", "lasso","elasticnet"), 3),
                    summary = rep(c("avg","pc"), each = 3),
                    model = rep(c("lasso","elasticnet","shim"), 2),
                    MoreArgs = list(x_train = result[["X_train"]],
                                    x_test = result[["X_test"]],
                                    y_train = result[["Y_train"]],
                                    y_test = result[["Y_test"]],
                                    stability = F,
                                    filter = F,
                                    filter_var = F,
                                    include_E = T,
                                    include_interaction = T,
                                    s0 = result[["S0"]],
                                    p = p,
                                    gene_groups = result[["clustersAll"]],
                                    clust_type = "clust"),
                    SIMPLIFY = F,
                    USE.NAMES = F)

# result %>% names
clust_res %>% unlist

includeStability = TRUE
if (includeStability) {
  clust_stab <- mapply(function(summary,
                                model) mapply(clust_fun,
                                              x_train = result[["X_train_folds"]],
                                              y_train = result[["Y_train_folds"]],
                                              MoreArgs = list(stability = T,
                                                              x_test = result[["X_test"]],
                                                              summary = summary,
                                                              model = model,
                                                              filter = F,
                                                              filter_var = F,
                                                              include_E = T,
                                                              include_interaction = T,
                                                              gene_groups = result[["clustersAll"]],
                                                              p = p,
                                                              clust_type = "clust"),
                                              SIMPLIFY = F),
                       #summary = rep(c("pc","spc","avg"), each = 3),
                       #model = rep(c("lm", "lasso","elasticnet"), 3),
                       summary = rep(c("avg","pc"), each = 3),
                       model = rep(c("lasso","elasticnet","shim"), 2),
                       SIMPLIFY = F,
                       USE.NAMES = F)
  
  
  # Make the combinations of list elements
  ll <- lapply(seq_along(clust_stab), function(i) combn(clust_stab[[i]], 2, simplify = F))
  
  
  clust_labels <- function(summary, model) {
    paste0("clust",paste0("_",summary),paste0("_",model),"_","yes_")
  }
  
  clust_labs <- mapply(clust_labels,
                       #summary = rep(c("pc","spc","avg"), each = 3),
                       #model = rep(c("lm", "lasso","elasticnet"), 3),
                       summary = rep(c("avg","pc"), each = 3),
                       model = rep(c("lasso","elasticnet","shim"), 2),
                       USE.NAMES = F)
  
  
  # Pairwise correlations of the model coefficients for each of the 10 CV folds
  clust_mean_stab <- lapply(seq_along(ll), function(j) {
    lapply(c("pearson","spearman"), function(i) {
      res <- mean(sapply(ll[[j]] , function(x) WGCNA::cor(x[[1]]$coef.est, x[[2]]$coef.est, method = i,
                                                          use = 'pairwise.complete.obs')), na.rm = TRUE)
      names(res) <- paste0(clust_labs[[j]], i)
      return(res)
    }
    )
  }
  )
  
  clust_mean_stab %>% unlist
  
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

print("done clust and regress interaction")


## ---- Ecluster-and-regress----

# we will treat the clusters as fixed i.e., even if we filter, or
# do cross validation, the group labels are predetermined by the
# above clustering procedure
# This method is based on clusters derived without accounting for the environment
# AND accoundting for the environment
# So we want to see if adding these clusters makes a difference from what people
# might usually do, i.e just based on correlation without E

message("starting Environment cluster and regress with interaction")
includeStability = TRUE
includeInteraction = TRUE

Eclust_res <- mapply(clust_fun,
                     #summary = rep(c("pc","spc","avg"), each = 3),
                     #model = rep(c("lm", "lasso","elasticnet"), 3),
                     summary = rep(c("avg","pc"), each = 3),
                     model = rep(c("lasso","elasticnet","shim"), 2),
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
                                     gene_groups = result[["clustersAddon"]],
                                     clust_type = "Eclust"),
                     SIMPLIFY = F,
                     USE.NAMES = F)

# Eclust_res %>% names
# options(digits = 2, scipen=999)
Eclust_res %>% unlist

if (includeStability) {
  Eclust_stab <- mapply(function(summary,
                                 model) mapply(clust_fun,
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
                                                               clust_type = "Eclust"),
                                               SIMPLIFY = F),
                        #summary = rep(c("pc","spc","avg"), each = 3),
                        #model = rep(c("lm", "lasso","elasticnet"), 3),
                        summary = rep(c("avg","pc"), each = 3),
                        model = rep(c("lasso","elasticnet","shim"), 2),
                        SIMPLIFY = F,
                        USE.NAMES = F)
  
  
  # Make the combinations of list elements
  ll <- lapply(seq_along(Eclust_stab), function(i) combn(Eclust_stab[[i]], 2, simplify = F))
  
  Eclust_labels <- function(summary, model) {
    paste0("Eclust",paste0("_",summary),paste0("_",model),"_","yes_")
  }
  
  Eclust_labs <- mapply(Eclust_labels,
                        #summary = rep(c("pc","spc","avg"), each = 3),
                        #model = rep(c("lm", "lasso","elasticnet"), 3),
                        summary = rep(c("avg","pc"), each = 3),
                        model = rep(c("lasso","elasticnet","shim"), 2),
                        USE.NAMES = F)
  
  
  # Pairwise correlations of the model coefficients for each of the 10 CV folds
  Eclust_mean_stab <- lapply(seq_along(ll), function(j) {
    lapply(c("pearson","spearman"), function(i) {
      res <- mean(sapply(ll[[j]] , function(x) WGCNA::cor(x[[1]]$coef.est, x[[2]]$coef.est, method = i,
                                                          use = "pairwise.complete.obs")), na.rm = TRUE)
      names(res) <- paste0(Eclust_labs[[j]], i)
      return(res)
    }
    )
  }
  )
  
  Eclust_mean_stab %>% unlist
  
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


## ---- penalization ----

print("starting penalization with interaction")

pen_res <- mapply(pen_fun,
                  #model = c("scad","mcp","lasso","elasticnet"),
                  model = c("lasso","elasticnet"),
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

pen_res %>% unlist()

if (includeStability) {
  pen_stab <- mapply(function(model) mapply(pen_fun,
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
                     model = c("lasso","elasticnet"),
                     SIMPLIFY = F,
                     USE.NAMES = F)
  
  
  # Make the combinations of list elements
  ll <- lapply(seq_along(pen_stab), function(i) combn(pen_stab[[i]], 2, simplify = F))
  
  
  pen_labels <- function(model) {
    paste0("pen_na_",model,"_yes_")
  }
  
  pen_labs <- mapply(pen_labels,
                     #model = c("scad","mcp","lasso","elasticnet","ridge"),
                     model = c("lasso","elasticnet"),
                     USE.NAMES = F)
  
  
  # Pairwise correlations of the model coefficients for each of the 10 CV folds
  pen_mean_stab <- lapply(seq_along(ll), function(j) {
    lapply(c("pearson","spearman"), function(i) {
      res <- mean(sapply(ll[[j]] , function(x) WGCNA::cor(x[[1]]$coef.est, x[[2]]$coef.est, method = i,
                                                          use = "pairwise.complete.obs")), na.rm = TRUE)
      names(res) <- paste0(pen_labs[[j]], i)
      return(res)
    }
    )
  }
  )
  
  pen_mean_stab %>% unlist
  
  # Jaccard index
  pen_jacc <- lapply(seq_along(ll), function(j) {
    res <- mean(sapply(ll[[j]] , function(x) {
      A = x[[1]][coef.est != 0]$Gene
      B = x[[2]][coef.est != 0]$Gene
      if (length(A)==0 | length(B)==0) 0 else length(intersect(A,B))/length(union(A,B))
    }), na.rm = TRUE)
    names(res) <- paste0(pen_labs[[j]],"jacc")
    return(res)
  })
  
  pen_jacc %>% unlist
}

print("done penalization with interaction")

## ---- not-used ----
pheatmap::pheatmap(TOM_all, color = viridis(100), clustering_method = "average",
                   show_rownames = F, show_colnames = F)
