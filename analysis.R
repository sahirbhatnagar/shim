##################################
# R source code file for analyzing NIHPD data
# Created by Sahir, May 24
# Updated:
# NOTE: This directory is more upto data than bouchard/scripts/ccareport
# I am now using version control in bouchard/cca/
# because it was getting too messy in ccareport
##################################

rm(list=ls())
source("functions.R")
source("packages.R")
source("data_cleaning.R")

brain_probes[1:10]
colnames(DT_with_pheno)[1:40]

DT_with_pheno[,"Subject_Gender",with=F] %>% str
DT_with_pheno[is.na(income)] %>% dim
DT_with_pheno[is.na(WASI_Full_Scale_IQ)] %>% dim
DT_with_pheno[is.na(Subject_Gender)] %>% dim
DT_with_pheno[is.na(Site_Location)] %>% dim


## ---- test-stat ----
registerDoMC(cores = 30)
res <- mclapply(brain_probes, function(i) {
  
  # i = brain_probes[1]
  dat <- DT_with_pheno[,c("WASI_Full_Scale_IQ","Subject_Gender","Site_Location",i),with=F]
  dat <- dat[!is.na(WASI_Full_Scale_IQ)]
  include_E <- FALSE
  include_interaction <- FALSE
  fit <- lm(as.formula(paste("WASI_Full_Scale_IQ ~",i, if (include_E) "+E", if (include_interaction) paste0("+",i,":E"),
                             "+ Subject_Gender + Site_Location")),
            data = dat)
  #fit %>% summary
  coef.index <- if (include_interaction) paste0(i,":E") else i
  
  data.frame(probe = i,
             pvalue = as.numeric(pt(abs(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
                                    df = fit$df.residual, lower.tail = F)*2),
             test.stat = as.numeric(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
             'mean' = mean(dat[,i, with = F][[1]]),
             'sd' = sd(dat[,i,with = F][[1]]), stringsAsFactors = FALSE)
}, mc.cores = 30)

uni_results <- do.call(rbind, res)
uni_results <- as.data.table(uni_results)


# uni_results[order(pvalue)[1:1000]]
# uni_results[, hist(sd)]
# uni_results[, hist(test.stat)]



## ---- filter-stats ----

with(uni_results,
     plot(rank(mean)/length(mean), -log10(pvalue), pch=16, cex=0.45))

par(mfrow = c(1,2))
with(uni_results,
     plot(rank(sd)/length(sd), -log10(pvalue), pch=16, cex=0.45))

trsf = function(n) log10(n+1)
plot(ecdf(trsf(uni_results$sd)), xlab=body(trsf), main="")
dev.off()


theta = seq(from=0, to=0.5, by=0.1)
pBH = filtered_p(filter=uni_results$mean, test=uni_results$pvalue, theta=theta, method="BH")
head(pBH)

rejection_plot(pBH, at="sample",
               xlim=c(0, 0.5), ylim=c(0, 15000),
               xlab="FDR cutoff (Benjamini & Hochberg adjusted p-value)", main="")

theta = seq(from=0, to=0.8, by=0.02)
rejBH = filtered_R(alpha = 0.1,
                   filter=uni_results$mean, 
                   test=uni_results$pvalue, 
                   theta=theta, method="BH")

plot(theta, rejBH, type="l",
     xlab=expression(theta), ylab="number of rejections")


theta = seq(from=0, to=0.99, by=0.02)
filterChoices = data.frame(
  `mean` = uni_results$mean,
  `geneID` = 1:nrow(uni_results),
  `sd` = uni_results$sd,
  `t.stat` = uni_results$test.stat
)
rejChoices = sapply(filterChoices, function(f)
  filtered_R(alpha=0.1, filter=f, test=uni_results$pvalue, theta=theta, method="BH"))

library("RColorBrewer")
myColours = brewer.pal(ncol(filterChoices), "Set1")
matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
        xlab=expression(theta), ylab="number of rejections")
legend("topleft", legend=colnames(filterChoices), fill=myColours)

rejChoices

theta = theta[which.max(rejChoices[,"mean"])]
pass = with(uni_results, mean > quantile(mean, theta))

h1 = hist(uni_results$pvalue[!pass], breaks=50, plot=FALSE)
h2 = hist(uni_results$pvalue[pass], breaks=50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


resFilt = uni_results[pass,]
orderInPlot = order(resFilt$pvalue)
showInPlot = (resFilt$pvalue[orderInPlot] <= 0.06)
alpha = 0.1

plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]),
     xlim=c(0,4000))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)

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




## ---- data ----

idAll <- as.character(DT_with_pheno[!is.na(WASI_Full_Scale_IQ)][["Subject_ID"]])
idExposed <- as.character(DT_with_pheno[!is.na(WASI_Full_Scale_IQ)][E==1][["Subject_ID"]])
idUnexposed <- as.character(DT_with_pheno[!is.na(WASI_Full_Scale_IQ)][E==0][["Subject_ID"]])

filtered_probes <- uni_results[order(test.stat,decreasing = T)[1:500]]$probe

# X <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ), 
#                    c("Subject_ID","age_binary","Subject_Gender","Site_Location",brain_probes), with = F]

X <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ), c("E",filtered_probes), with = F]
dim(X)
# dimnames(X)[[1]] <- idAll
head(X)
str(X)

X0 <- X[E==0]
dim(X0)
# dimnames(X0)[[1]] <- idUnexposed

X1 <- X[E==1]
dim(X1)
# dimnames(X1)[[1]] <- idExposed


TOM_all <- TOMsimilarityFromExpr(X[,filtered_probes,with=F], nThreads = 35)
TOM_exposed <- TOMsimilarityFromExpr(X1[,filtered_probes,with=F], nThreads = 35)
TOM_unexposed <- TOMsimilarityFromExpr(X0[,filtered_probes,with=F], nThreads = 35)
TOM_diff <- abs(TOM_exposed-TOM_unexposed)

Y <- as.matrix(DT_with_pheno[!is.na(WASI_Full_Scale_IQ),"WASI_Full_Scale_IQ",with=F])
dimnames(Y)[[1]] <- idAll


## ---- data-prep-for-analysis ----

# partition-data
trainIndex <- caret::createDataPartition(Y[,"WASI_Full_Scale_IQ"], p = 2/3, 
                                         list = FALSE, times = 1)

# DT_train <- DT_with_pheno[trainIndex[,1],,]
# DT_test <- DT[-trainIndex,]

# X_train and X_test contain the environment variable
X_train <- X[trainIndex[,1]]
dim(X_train)
Y_train <- Y[trainIndex[,1],]
str(Y_train)

# fit.cv <- cv.glmnet(as.matrix(X_train), Y_train)
# plot(fit.cv)

X_test <- X[-trainIndex[,1]]
dim(X_test)
Y_test <- Y[-trainIndex[,1]]

mse_null <- crossprod(mean(Y_test[[1]]) - Y_test[[1]])/length(Y_test[[1]])

# training gene expression data, this does not include the environment variable
genes_e0 <- as.matrix(X_train[E==0,filtered_probes,with=F])
dim(genes_e0)
genes_e1 <- as.matrix(X_train[E==1,filtered_probes,with=F])
dim(genes_e1)
genes_all <- rbind(genes_e0,genes_e1)

message("Calculating similarity matrices")

# test set gene expression data
genes_all_test <- as.matrix(X_test[,filtered_probes, with = F])
dim(genes_all_test)

# corr_train_e0 <- WGCNA::cor(genes_e0)
# corr_train_e1 <- WGCNA::cor(genes_e1)
# corr_train_diff <- abs(corr_train_e1 - corr_train_e0)
corr_train_all <- WGCNA::cor(genes_all)

tom_train_e0 <- WGCNA::TOMsimilarityFromExpr(genes_e0)
dim(tom_train_e0)
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

message("Creating CV folds from training data")

# Folds for Cross validation
folds_train <- caret::createFolds(Y_train, k = 10, list = T)
# DT_train_folds <- lapply(folds_train, function(i) DT_train[-i,])
# these have the environment variable
X_train_folds <- lapply(folds_train, function(i) X_train[-i])
Y_train_folds <- lapply(folds_train, function(i) Y_train[-i])

cluster_distance <- "tom"
clustMethod <- "hclust"
cutMethod <- "dynamic"
method <- "complete"
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

## ---- not-used ----
pheatmap::pheatmap(TOM_all, color = viridis(100), clustering_method = "average",
                   show_rownames = F, show_colnames = F)
