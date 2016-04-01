##################################
# R source code file for simulating data, WGCNA style
#
# Created by Sahir,  March 25, 2016
# Updated:
# Notes:
# This is based on code from Network analysis book by Horvath
##################################

library(viridis)
library(pheatmap)
library(WGCNA)
library(magrittr)
library(doMC)
library(Matrix)
library(data.table)
library(protoclust)

dev.off()
rm(list = ls())

## ---- functions ----

#' @param rho correlations between green module and temp vector T, and red
#'   module and temp vector T. A numeric vector of length 2
simModule <- function(n, p, rho, exposed, ...) {
  
  if (exposed) {
    #Step 1: simulate the seed module eigengenes
    sMEturquoise <- rnorm(n)
    
    #expected cor(sMEblue,sMEturquoise) = 0.60
    sMEblue <- 0.60 * sMEturquoise + sqrt(1 - 0.60 ^ 2) * rnorm(n)
    
    sMEyellow <- rnorm(n)
    
    sMEgreen <- rnorm(n)
    
    #expected cor(e.continuous,seed.ME)=0.95
    temp0 <- rho[1] * sMEgreen + sqrt(1 - rho[1] ^ 2) * rnorm(n)
    
    #expected cor(y.continuous,seed.ME) <- -0.95
    sMEred <- rho[2] * temp0 + sqrt(1 - rho[2] ^ 2) * rnorm(n)
    
    datsME <- data.frame(sMEturquoise,sMEblue,sMEred,sMEgreen,sMEyellow) 
    
    dat1 <- WGCNA::simulateDatExpr(eigengenes = datsME, nGenes = p, ...)
  } else {
    
    #Step 1: simulate the seed module eigengenes
    sMEturquoise <- rnorm(n)
    
    #expected cor(sMEblue,sMEturquoise) = 0.60
    sMEblue <- 0.60 * sMEturquoise + sqrt(1 - 0.60 ^ 2) * rnorm(n)
    
    sMEyellow <- rnorm(n)
    
    sMEgreen <- rnorm(n)
    
    #expected cor(e.continuous,seed.ME)=0.95
    temp0 <- rho[1] * sMEgreen + sqrt(1 - rho[1] ^ 2) * rnorm(n)
    
    #expected cor(y.continuous,seed.ME) <- -0.95
    sMEred <- rho[2] * temp0 + sqrt(1 - rho[2] ^ 2) * rnorm(n)
    
    datsME <- data.frame(sMEturquoise,sMEblue,sMEred,sMEgreen,sMEyellow)
    
    dat1 <- WGCNA::simulateDatExpr(eigengenes = datsME, nGenes = p, ...)
    
  } 
  
  return(dat1)
}


#' Function to generate heatmap
#' @decription x matrix of true correlation (P x P matrix where P is the number 
#'   of genes). Must be object of class similarity

plot.similarity <- function(x,
                            color = viridis(100),
                            truemodule, ...){
  
  annotation_col <- data.frame(
    module = factor(truemodule, 
                    labels = c("Grey","Turquoise","Blue","Red",
                               "Green","Yellow")))
  
  rownames(annotation_col) <- dimnames(x)[[2]]
  ann_colors <- list(
    module = c(Turquoise = "turquoise",
               Blue = "blue",
               Red = "red",
               Green = "green",
               Yellow = "yellow",
               Grey = "grey90")
  )
  
  pheatmap(x,
           show_rownames = F, show_colnames = F,
           color = color,
           annotation_col = annotation_col,
           annotation_row = annotation_col,
           annotation_colors = ann_colors,
           annotation_names_row = FALSE,
           annotation_names_col = FALSE,
           drop_levels = FALSE,
           annotation_legend = FALSE, ...)
  
  
  # breaks = seq(min(min_max_heat$V2), max(min_max_heat$V1), length.out = 101) ,
  # legend_breaks = round(seq(min(min_max_heat$V2), max(min_max_heat$V1), 
  # length.out = 12),1),
  # legend_labels = round(seq(min(min_max_heat$V2), max(min_max_heat$V1), 
  # length.out = 12),1),
  # drop_levels = FALSE, ...)
}


proto1 <- function(x,k) list(cluster = {
  as.numeric(protoclust::protocut(protoclust(as.dist(x)), k = k)$cl)})

hclustWard <- function(x,k) list(cluster = {
  as.numeric(cutree(hclust(as.dist(x), method = "ward.D2"), k = k))
})


bigcorPar <- function(data.all, data.e0, data.e1, alpha = 1.5, threshold = 1, 
                      nblocks = 10, ncore = 2, ...){
  
  # data.all = X
  # data.e0 = X0
  # data.e1 = X1
  # nblocks = 10
  # alpha = 1.5; threshold = 1
  
  NCOL <- ncol(data.all)
  
  SPLIT <- split(1:NCOL, ceiling(seq_along(1:NCOL)/nblocks))
  
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  result <- foreach(i = 1:nrow(COMBS), .combine = rbind) %dopar% {
    #i=2567
    #COMBS[1:3000,]
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    
    rho.all <- tril(WGCNA::corFast(data.all[, G1], data.all[, G2]), k = -1)
    rho.gd <- tril(WGCNA::corFast(data.e1[, G1], data.e1[, G2]),  k = -1)
    rho.ngd <- tril(WGCNA::corFast(data.e0[, G1], data.e0[, G2]), k = -1)
    
    S <- abs(rho.gd + rho.ngd - alpha * rho.all)
    arr <- which(S > threshold, arr.ind = T)
    as.data.table(
      data.frame(
        cbind(
          "score" = S[arr] ,
          "gene1" = dimnames(rho.all[arr[,1],arr[,2], drop = F])[[1]],
          "gene2" = dimnames(rho.all[arr[,1],arr[,2], drop = F])[[2]]),
        stringsAsFactors = FALSE))
  }
  
  return(result)
}

## ---- data ----
n0 = 50
n1 = 50
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

dim(X)

corrX <- cor(X)

X0 <- X[paste0("Subject",1:n0),]
X1 <- X[paste0("Subject",(n0+1):n),]

class(corrX) <- c("similarity",class(corrX))

corrX0 <- cor(X0)
class(corrX0) <- c("similarity",class(corrX0))

corrX1 <- cor(X1)
class(corrX1) <- c("similarity",class(corrX1))

diffCorr <- abs(corrX1-corrX0)
class(diffCorr) <- c("similarity",class(diffCorr))

TOMX <- TOMsimilarityFromExpr(X)
class(TOMX) <- c("similarity",class(TOMX))

TOMX0 <- TOMsimilarityFromExpr(X0)
class(TOMX0) <- c("similarity",class(TOMX0))

TOMX1 <- TOMsimilarityFromExpr(X1)
class(TOMX1) <- c("similarity",class(TOMX1))

diffTOM <- abs(TOMX1 - TOMX0)
class(diffTOM) <- c("similarity",class(diffTOM))

alpha <- 1.5
Scorr <- abs(corrX0 + corrX1 - alpha * corrX)
class(Scorr) <- c("similarity", class(Scorr))

Stom <- abs(TOMX0 + TOMX1 - alpha * TOMX)
class(Stom) <- c("similarity", class(Stom))

## ---- heat-corr-all ----
hc <- hclust(as.dist(1 - corrX), method = "ward.D2")
plot(corrX, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-corr-e0 ----
hc <- hclust(as.dist(1 - corrX0), method = "ward.D2")
plot(corrX0, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-corr-e1 ----
hc <- hclust(as.dist(1 - corrX1), method = "ward.D2")
plot(corrX1, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-corr-diff ----
hc <- hclust(as.dist(diffCorr), method = "ward.D2")
plot(diffCorr, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-tom-all ----
hc <- hclust(as.dist(dissTOMX), method = "ward.D2")
plot(TOMsimilarityFromExpr(X), truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-tom-e0 ----
plot(dissTOMX0, truemodule = truemodule1, cluster_rows = F, cluster_cols = F)

## ---- heat-tom-e1 ----
plot(dissTOMX1, truemodule = truemodule1, cluster_rows = F, cluster_cols = F)

## ---- heat-tom-diff ----
plot(diffTOM, truemodule = truemodule1, cluster_rows = T, cluster_cols = T)

## ---- cor-scor ----
plot(Scorr, truemodule = truemodule1, cluster_rows = T, cluster_cols = T)

## ---- cor-scor-tom ----
plot(Stom, truemodule = truemodule1, cluster_rows = T, cluster_cols = T)

## ---- cluster-cor-scor ----

## ---- module-eigengenes ----

MEList <- blockwiseModules(X)

## ---- gap-statistic ----




gapScorr <- cluster::clusGap(dissTOMX, FUNcluster = proto1, K.max = 20, B = 50)
gapScorr
plot(gapScorr, main = "clusGap(., FUN = protoclust, n.start=20, B= 60)")

lapply(c("firstSEmax", "Tibs2001SEmax", "globalSEmax",
         "firstmax", "globalmax"), function(i) maxSE(f = gapScorr$Tab[, "gap"], 
                                                     SE.f = gapScorr$Tab[, "SE.sim"], 
                                                     method = i, 
                                                     SE.factor = 1)) 


hc <- protoclust(as.dist(1 - corrX))

hc <- protoclust(as.dist(dissTOMX))
class(hc) <- "hclust"


dissTOMX %>% str
dimnames(dissTOMX)[[1]] <- paste0("Gene",1:500)
dimnames(dissTOMX)[[2]] <- paste0("Gene",1:500)

plot(dissTOMX, truemodule = truemodule1, 
     cluster_rows = hc, 
     cluster_cols = hc,
     cutree_cols = 3
)

# ======================================#

hc <- protoclust(as.dist(diffCorr))
dev.off()
class(hc) <- "hclust"
as.dist(diffCorr) %>% str

hc <- hclust(as.dist(diffCorr), method = "ward.D2")
plot(hc)
labelDynamicHybrid <- cutreeDynamic(hc, distM = diffCorr,
                                    cutHeight = 0.995, deepSplit = 3,
                                    pamRespectsDendro = FALSE,
                                    minClusterSize = 20)

labelDynamicHybrid %>% unique %>% length

plot(Scorr, truemodule = truemodule1,cluster_rows = hc, cluster_cols = hc,
     cutree_cols = 8)
plot(Scorr, truemodule = truemodule1,cluster_rows = T, cluster_cols = T,
     clustering_callback = callback, cutree_cols = 5)
# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap()

gapDiff <- cluster::clusGap(1-corrX, FUNcluster = hclustWard, K.max = 20, B = 20)
print(gapDiff, method = "Tibs")
gapDiff

maxSE(f = gapDiff$Tab[, "gap"], SE.f = gapDiff$Tab[, "SE.sim"], 
      method = "Tibs2001SEmax", 
      SE.factor = 1)

plot(gapDiff, main = "clusGap(., FUN = kmeans, n.start=20, B= 60)")
hc <- hclust(as.dist(1 - corrX), method = "ward.D2")
plot(hc)

plot(corrX, truemodule = truemodule1, 
     cluster_rows = hc, 
     cluster_cols = hc,
     cutree_cols = 5
)




## ---- proto-clust ----

# generate some data:
dev.off()
set.seed(1)
n <- 100
p <- 2
x <- matrix(rnorm(n * p), n, p)
rownames(x) <- paste("A", 1:n, sep="")
d <- dist(x)
d <- S
# perform minimax linkage clustering:
hc <- protoclust(d)
# cut the tree to yield a 10-cluster clustering:
k <- 10 # number of clusters
cut <- protocut(hc, k=k)
h <- hc$height[n - k]
# plot dendrogram (and show cut):
plotwithprototypes(hc, imerge=cut$imerge, col=2)
abline(h=h, lty=2)

## ---- module-preserved ----

# We now set up the multi-set expression data
# and corresponding module colors:
setLabels = c("E0", "E1")
multiExpr=list(E0=list(data=X0),
               E1=list(data=X1))
multiColor=list(E0 = moduleColors, E1=moduleColors)
nPermutations1=3
mp = modulePreservation(multiExpr, multiColor,
                        referenceNetworks = 1, nPermutations = nPermutations1,
                        randomSeed = 1, quickCor = 0, verbose = 3)
# specify the reference and the test networks
ref=1; test = 2
Obs.PreservationStats= mp$preservation$observed[[ref]][[test]]
Z.PreservationStats=mp$preservation$Z[[ref]][[test]]
# Look at the observed preservation statistics
Obs.PreservationStats
# Z statistics from the permutation test analysis
Z.PreservationStats
# Let us now visualize the data.
modColors = rownames(Obs.PreservationStats)
moduleSize = Obs.PreservationStats$moduleSize
# we will omit the grey module (background genes)
# and the gold module (random sample of genes)
selectModules = !(modColors %in% c("grey", "gold"))
# Text labels for points
point.label = modColors[selectModules]
#Composite preservation statistics
medianRank=Obs.PreservationStats$medianRank.pres
Zsummary=Z.PreservationStats$Zsummary.pres
par(mfrow=c(2,2),mar = c(4.5,4.5,2.5,1))
# plot medianRank versus module size
plot(moduleSize[selectModules],medianRank[selectModules],col=1,
     bg=modColors[selectModules],pch = 21,main="medianRank
     Preservation",
     cex = 2, ylab ="medianRank",xlab="Module size", log="x")
labelPoints(moduleSize[selectModules],medianRank[selectModules],
            point.label,cex=1,offs=0.03)
# plot Zsummary versus module size
plot(moduleSize[selectModules],Zsummary[selectModules], col = 1,
     bg=modColors[selectModules],pch = 21,
     main="Zsummary Preservation",
     cex=2,ylab ="Zsummary", xlab = "Module size", log = "x")
labelPoints(moduleSize[selectModules],Zsummary[selectModules],
            point.label,cex=1,offs=0.03)
# Add threshold lines for Zsummary
abline(h=0); abline(h=2, col = "blue", lty = 2)
abline(h=10, col = "red", lty = 2)

## ---- not-used ----


# We now calculate the weighted adjacency matrix, using the
# power 6
A = adjacency(X1, power = 6)
#define a dissimilarity based on the topological overlap
dissTOM =TOMdist(A)
#hierarchical clustering
geneTree = hclust(as.dist(dissTOM),method="average")

plot(geneTree)

# here we define the modules by cutting branches
moduleLabelsManual1 = cutreeDynamic(dendro=geneTree,distM=dissTOM,
                                    method="hybrid",deepSplit=2,pamRespectsDendro=F,
                                    minClusterSize=30)

# Convert labels to colors for plotting
moduleColorsManual2=labels2colors(moduleLabelsManual1)



# Set the diagonal of the TOM disscimilarity to NA
diag(dissTOM) = NA
# Transform dissTOM with a power to enhance visibility
TOMplot(dissTOM ^ 7,geneTree,moduleColorsManual2,
        main = "Network heatmap plot, all genes")



dim(datExpr)

# simulated expression data
#datExpr <- dat1$datExpr
datExpr <- rbind(d0$datExpr, d1$datExpr)
# next we assign names to the rows and columns of datExpr
datExpr <- data.frame(datExpr)
ArrayName <- paste("Sample",1:dim(datExpr)[[1]],sep="")
dimnames(datExpr)[[1]] <- ArrayName
GeneName <- paste("Gene",1:dim(datExpr)[[2]],sep="")
dimnames(datExpr)[[2]] <- GeneName

#simulated true module color of each gene
truemodule <- dat1$setLabels

table(truemodule)

## ---- simulate-data ----

#set the seed of the random number generator
#set.seed(1)

# number of samples
n <- 100

# number of genes
p <- 1000

#Step 1: simulate the seed module eigengenes
sMEturquoise <- rnorm(n)

#expected cor(sMEblue,sMEturquoise) = 0.75
sMEblue <- 0.75 * sMEturquoise + sqrt(1 - 0.75 ^ 2) * rnorm(n)

cor(sMEblue, sMEturquoise)

sMEyellow <- rnorm(n)
sMEgreen <- rnorm(n)

#expected cor(e.continuous,seed.ME)=0.95
e.continuous <- 0.75 * sMEgreen + sqrt(1 - 0.75 ^ 2) * rnorm(n)

#expected cor(y.continuous,seed.ME) <- -0.95
sMEred <- 0.75 * e.continuous + sqrt(1 - 0.75 ^ 2) * rnorm(n)

cor(e.continuous, sMEgreen)
cor(e.continuous, sMEred)
cor(sMEred, sMEgreen)

#dichotomize y.continuous
e <- ifelse(e.continuous > median(e.continuous),2,1)
table(e)

datsME <- data.frame(e,sMEturquoise,sMEblue,sMEred,sMEgreen,sMEyellow)
signif(cor(datsME,use = "p"),1)
datsME <- data.frame(sMEturquoise,sMEblue,sMEred,sMEgreen,sMEyellow)

## ---- weighted-coex-network ----

# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1 <- abs(cor(datExpr,use="p"))^6
#ADJ1 <- WGCNA::TOMsimilarityFromExpr(datExpr)
# When you have relatively few genes (<5000) use the following code
k <- as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k <- softConnectivity(datE=datExpr,power=6)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow = c(1,2))
hist(k)
scaleFreePlot(k, main = "Check scale free topology\n")

datsME %>% colnames

annotation_col <- data.frame(
  module = factor(truemodule, labels = c("Grey","Turquoise","Blue","Red","Green","Yellow")))

rownames(annotation_col) <- dimnames(datExpr)[[2]]
ann_colors <- list(
  module = c(Turquoise = "turquoise",
             Blue = "blue",
             Red = "red",
             Green = "green",
             Yellow = "yellow",
             Grey = "grey90")
)

pheatmap(cor(datExpr,use = "p"),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F, show_colnames = F,
         color = viridis(100),
         annotation_col = annotation_col,
         annotation_row = annotation_col,
         annotation_colors = ann_colors)

library(cluster)
#Turn the similarity measure into a dissimilarity measure
dissTOM <- TOMdist(ADJ1)
# define a hierarchical tree with hclust or with flashClust
hierTOM <- hclust(as.dist(dissTOM), method = "average")
plot(hierTOM)

# Here we use static branch cutting with height 0.995
# and minimum cluster size of 30.
colorStaticTOM <- cutreeStaticColor(hierTOM, cutHeight = 0.995,
                                    minSize = 30)

# We now use two Dynamic Tree Cut methods in which the height
# cut-off plays a minor role. The first method is called the tree
# method and only uses the dendrogram as input.
branch.number <- cutreeDynamic(hierTOM, method = "tree", cutHeight = 0.995,
                               deepSplit = F, minClusterSize = 30)

# This function transforms the branch numbers into colors
colorDynamicTOM <- labels2colors(branch.number)

labelDynamicHybrid <- cutreeDynamic(hierTOM, distM = dissTOM,
                                    cutHeight = 0.995, deepSplit = 1,
                                    pamRespectsDendro = FALSE,
                                    minClusterSize = 30)

colorDynamicHybridTOM <- labels2colors(labelDynamicHybrid)

# Plot results of all module detection methods together:
plotDendroAndColors(dendro = hierTOM,
                    colors = data.frame(truemodule, colorStaticTOM, 
                                        colorDynamicTOM,
                                        colorDynamicHybridTOM),
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main ="Cluster dendrogram and color-band")

randIndex(table(colorStaticTOM,truemodule),adjust=F)

randIndex(table(colorDynamicTOM,truemodule),adjust=F)

randIndex(table(colorDynamicHybridTOM,truemodule),adjust=F)