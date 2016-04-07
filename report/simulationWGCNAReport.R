##################################
# R source code file for simulating data, WGCNA style. This
# is only used for the report. See funcions.R and data.R for actual simulation
# source code
# Git: this is on the eclust repo, simulation branch
# Created by Sahir,  March 25, 2016
# Updated: April 2, 2016
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
library(cluster)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(factoextra)


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

fisherTransform <- function (n1, r1, n2, r2) {
  num1a <- which(r1 >= 0.99)
  num2a <- which(r2 >= 0.99)
  r1[num1a] <- 0.99
  r2[num2a] <- 0.99
  num1b <- which(r1 <= -0.99)
  num2b <- which(r2 <= -0.99)
  r1[num1b] <- -0.99
  r2[num2b] <- -0.99
  # atanh (inverse hyperbolic tangent) simplifies to
  # 0.5 * log(1+r)/log(1-r) , for r < 1
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  dz <- (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
  pv <- 2 * (1 - pnorm(abs(dz)))
  return(list(diff = dz, pval = pv))
}

#' Calculate Fisher's Z test for correlations
fisherZ <- function(n0, cor0, n1, cor1) {

  # n0 = 50
  # n1 = 50
  # cor0 = corrX0
  # cor1 = corrX1

  # by default this doesnt include the diagonal
  # this collapses the correlation matrix by columns
  ccc0 <- as.vector(cor0[lower.tri(cor0)])
  ccc1 <- as.vector(cor1[lower.tri(cor1)])

  p <- nrow(cor1)

  # number of Z statistics to calculate (p choose 2)
  geneNames <- rownames(cor1)

  zstat <- fisherTransform(n0, ccc0, n1, ccc1)$diff

  # convert vector to symmetric matrix
  zMat <- diag(p)
  zMat[lower.tri(zMat)] <- zstat
  zMat <- zMat + t(zMat) - diag(diag(zMat))
  dimnames(zMat) <- list(geneNames,geneNames)
  class(zMat) <- c("similarity", class(zMat))
  return(zMat)
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
  as.numeric(cutree(hclust(as.dist(x), method = "complete"), k = k))
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


#' Modified version of WGCNA::moduleEigengenes
#'
#' @description this function retunrs the PC instead of just the loading
#' as well as the mean and sd for each module, that will be used in predicting
#' the PCs for the test data
firstPC <- function (expr, colors, impute = TRUE, nPC = 1, align = "along average",
          excludeGrey = FALSE, grey = if (is.numeric(colors)) 0 else "grey",
          subHubs = TRUE, trapErrors = FALSE, returnValidOnly = trapErrors,
          softPower = 6, scale = TRUE, verbose = 0, indent = 0) {

  # expr = expr; colors = clusters$module; impute = TRUE; nPC = 1; align = "along average";
  # excludeGrey = FALSE; grey = if (is.numeric(colors)) 0 else "grey";
  # subHubs = TRUE; trapErrors = FALSE; returnValidOnly = trapErrors;
  # softPower = 6; scale = TRUE; verbose = 0; indent = 0;

  if (is.null(expr)) {
    stop("moduleEigengenes: Error: expr is NULL. ")
  }
  if (is.null(colors)) {
    stop("moduleEigengenes: Error: colors is NULL. ")
  }
  if (is.null(dim(expr)) || length(dim(expr)) != 2)
    stop("moduleEigengenes: Error: expr must be two-dimensional.")
  if (dim(expr)[2] != length(colors))
    stop("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene).")
  if (is.factor(colors)) {
    nl = nlevels(colors)
    nlDrop = nlevels(colors[, drop = TRUE])
    if (nl > nlDrop)
      stop(paste("Argument 'colors' contains unused levels (empty modules). ",
                 "Use colors[, drop=TRUE] to get rid of them."))
  }

  # maxVarExplained = 10
  # if (nPC > maxVarExplained)
  #   warning(paste("Given nPC is too large. Will use value",
  #                 maxVarExplained))
  #
  # nVarExplained = min(nPC, maxVarExplained)

  modlevels = levels(factor(colors))

  if (excludeGrey)
    if (sum(as.character(modlevels) != as.character(grey)) >
        0) {
      modlevels = modlevels[as.character(modlevels) !=
                              as.character(grey)]
    } else {
    stop(paste("Color levels are empty. Possible reason: the only color is grey",
               "and grey module is excluded from the calculation."))
  }

  # these are the loadings aka the first eigenvector for each module
  # length of these vectors will vary depending on the size of the module
  eigenVectors <- vector("list", length(modlevels))

  # these are the actual PC's aka the data %*% eigenvector
  # each column will be a n-dimensional vector.. i.e. a value for each person
  PC <- data.frame(matrix(NA, nrow = dim(expr)[[1]],
                                ncol = length(modlevels)))

  # list to store prcomp objects
  prcompObj <- vector("list", length(modlevels))


  # this is the average expression in a module for each subject
  # so this is a n x length(modlevels) matrix
  averExpr <- data.frame(matrix(NA, nrow = dim(expr)[[1]],
                               ncol = length(modlevels)))

  varExpl <- vector("double", length(modlevels))

  # validMEs = rep(TRUE, length(modlevels))
  # validAEs = rep(FALSE, length(modlevels))

  # these are the means and sds used for subsequent predictions
  means = vector("list", length(modlevels))
  sds = vector("list", length(modlevels))

  # isPC = rep(TRUE, length(modlevels))
  # isHub = rep(FALSE, length(modlevels))
  validColors = colors

  # names(eigenVectors) = paste(moduleColor.getMEprefix(), modlevels,
  #                          sep = "")
  names(PC) = paste("PC", modlevels, sep = ":")
  names(averExpr) = paste("AE", modlevels, sep = ":")

  for (i in c(1:length(modlevels))) {
    #i=1
    if (verbose > 1)
      printFlush(paste(spaces, "moduleEigengenes : Working on ME for module",
                       modlevels[i]))
    modulename = modlevels[i]
    restrict1 = as.character(colors) == as.character(modulename)
    if (verbose > 2)
      printFlush(paste(spaces, " ...", sum(restrict1),
                       "genes"))

    datModule <- as.matrix(expr[, restrict1])

    # dim(datModule)
    # dim(t(datModule))
    # dim(expr)

    # using prcomp first (need to use untransposed data!)
    prcompObj[[i]] <- prcomp(datModule, center = scale, scale. = scale)
    #View(stats:::prcomp.default)
    # prcompObj[[i]]$x %>% dim
    # prcompObj[[i]] %>% names
    # prcompObj[[i]]$rotation %>% dim

    eigenVectors[[i]] <- prcompObj[[i]]$rotation[,1, drop = F]

    averExpr[,i] <- rowMeans(scale(datModule, center = scale, scale = scale),
                             na.rm = TRUE)

    varExpl[[i]] <- factoextra::get_eigenvalue(prcompObj[[i]])[1,"variance.percent"]
    # corAve = cor(averExpr[,i], prcompObj[[i]]$rotation[,1],
    #              use = "p")
    # if (!is.finite(corAve)) corAve = 0
    # if (corAve < 0) prcompObj[[i]]$rotation[,1] = -prcompObj[[i]]$rotation[,1]

    PC[, i] <- predict(prcompObj[[i]])[,1]
    # plot(PC[, i], prcompObj[[i]]$x[,1])
    #means[i] <- prcompObj[[i]]$center
    #sds[i] <- prcompObj[[i]]$scale


  }

  list(eigengenes = eigenVectors, averageExpr = averExpr,
       varExplained = varExpl, validColors = validColors,
       PC = PC, prcompObj = prcompObj)
}



#' Cluster similarity matrix and return cluster membership of each gene
#'
#' @param x similarity matrix. must have non-NULL dimnames i.e., the rows and
#'   columns should be labelled preferable "Gene1, Gene2, ..."
#' @param distanceMethod  one of "euclidean","maximum","manhattan", "canberra",
#'   "binary","minkowski" to be passed to \code{dist} function. If missing, then
#'   this function will take 1-x as the dissimilarity measure. This
#'   functionality is for diffCorr and diffTOM matrices which need to be
#'   converted to a distance type matrix.
#' @param clusterMethod how to cluster the data
#' @param cutMethod what method to use to cut the dendrogram. \code{dynamic}
#'   refers to dynamicTreeCut library, \code{gap} is Tibshirani's gap statistic
#'   \code{fixed} is a fixed number specified by the \code{nClusters} argument
#' @param nClusters number of clusters. Only used if \code{cutMethod = fixed}
#' @param method the agglomeration method to be used. This should be (an
#'   unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'   "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC)
#'   or "centroid" (= UPGMC).
clusterSimilarity <- function(x,
                              expr,
                              distanceMethod,
                              clustMethod = c("hclust", "protoclust"),
                              cutMethod = c("dynamic","gap", "fixed"),
                              nClusters,
                              method = c("complete", "average", "ward.D2",
                                         "single", "ward.D", "mcquitty",
                                         "median", "centroid"),
                              summary = c("avg", "pca"),
                              K.max = 10, B = 50) {

  x = corrX ; expr = X
  dim(X)
  clustMethod = c("hclust")
  cutMethod = c("dynamic")
  nClusters = 6
  method = c("complete")
  summary = c("pca")
  K.max = 10; B = 50
  distance = as.dist(1 - x)

  geneNames <- dimnames(x)[[1]]
  p <- nrow(x)
  method <- match.arg(method)
  cutMethod <- match.arg(cutMethod)
  clustMethod <- match.arg(clustMethod)

  distance <- if (missing(distanceMethod)) {
    as.dist(1 - x)
  } else dist(x = x, method = distanceMethod)

  hc <- switch(clustMethod,
               hclust = {
                 hclust(distance, method = method)
               },
               protoclust = {
                 protoclust(distance)
               }
  )

  #plot(hc)
  # create cluster function used if Gap statistic is requested
  # its as.dist(x) here because I am passing the
  # 1-x matrix to the cluster::clusGap function
  if (cutMethod == "gap") {

    FUNcluster <- if (missing(distanceMethod)) {
      switch(clustMethod,
             hclust = {
               function(xMat,k) list(cluster = {
                 as.numeric(
                   cutree(
                     hclust(as.dist(xMat), method = method), k = k
                   )
                 )
               })
             },
             protoclust = {
               function(xMat,k) list(cluster = {
                 as.numeric(protoclust::protocut(protoclust(as.dist(xMat)),
                                                 k = k)$cl)})
             }
      )
    } else {
      switch(clustMethod,
             hclust = {
               function(xMat,k) list(cluster = {
                 as.numeric(cutree(hclust(dist(xMat, method = distanceMethod),
                                          method = method), k = k))})
             },
             protoclust = {
               function(xMat,k) list(cluster = {
                 as.numeric(protoclust::protocut(
                   protoclust(dist(xMat, method = distanceMethod)),
                   k = k)$cl)})
             }
      )
    }
    #return(FUNcluster)
  }


  clustAssignment <- switch(cutMethod,
                            dynamic = {
                              if (clustMethod == "hclust") {
                                dynamicTreeCut::cutreeDynamic(
                                  hc,
                                  distM = as.matrix(distance),
                                  cutHeight = 0.995, deepSplit = 1,
                                  pamRespectsDendro = FALSE,
                                  minClusterSize = 20)
                              } else {
                                hcMod <- hc
                                class(hcMod) <- "hclust"
                                dynamicTreeCut::cutreeDynamic(
                                  hcMod,
                                  distM = as.matrix(distance),
                                  cutHeight = 0.995, deepSplit = 1,
                                  pamRespectsDendro = FALSE,
                                  minClusterSize = 20)
                              }
                            },
                            gap = {
                              if (clustMethod == "hclust") {
                                gapResult <- cluster::clusGap(1 - x,
                                                              FUNcluster = FUNcluster,
                                                              K.max = K.max,
                                                              B = B)
                                nClustGap <- maxSE(f = gapResult$Tab[, "gap"],
                                                   SE.f = gapResult$Tab[, "SE.sim"],
                                                   method = "Tibs2001SEmax",
                                                   SE.factor = 1)
                                cutree(hc, nClustGap)

                              } else {
                                gapResult <- cluster::clusGap(1 - x,
                                                              FUNcluster = FUNcluster,
                                                              K.max = K.max,
                                                              B = B)
                                nClustGap <- maxSE(f = gapResult$Tab[, "gap"],
                                                   SE.f = gapResult$Tab[, "SE.sim"],
                                                   method = "Tibs2001SEmax",
                                                   SE.factor = 1)
                                protocut(hc, k = nClustGap)[["cl"]]
                              }
                            },
                            fixed = {
                              if (clustMethod == "hclust") {
                                cutree(hc, nClusters)
                              } else protocut(hc, k = nClusters)[["cl"]]
                            }
  )

  # check if all cluster groups are 0 which means no cluster
  # assignment and everyone is in their own group
  plot(clustAssignment)
  clusters <- data.table(gene = geneNames,
                         cluster = if (all(clustAssignment == 0))
                           1:p else clustAssignment)
  #setkey(clusters, "cluster")

  # convert cluster numbers to colors which define modules
  clusters[, module := WGCNA::labels2colors(cluster)]
  clusters[, table(cluster,module)]


  # note that the align argument acts as follows if equal to "along average"
  # which is the default: it take the correlation between the average expression
  # in a module and the 1st eigenvector in a module and checks if its less
  # than 0, if its less than 0, then the moduleEigengenes function multiplies
  # the 1st eigenvector by -1, else it returns the unmodified 1st eigenvector
  # note that moduleEigengenes function returns the 1st eigenvector which is
  # equivalent to the rotation returned by prcomp, and what is used in
  # predict.prcomp to calculate the actual PCs.
  # to calculate PC's the following are all equivalent:
  # all.equal((expr %*% prcomp.object$rotation)[,1],
  # predict(prcomp.object)[,1],prcomp.object$x[,1])
  #
  # these are equivalent
  # p <- WGCNA::moduleEigengenes(expr = expr[, clusters$gene],
  #                              colors = clusters$module,
  #                              align = "",
  #                              scale = FALSE)
  # l <- prcomp(t(expr[, which(clusters$module %in% "blue")]), scale. = FALSE,
  # center = FALSE)
  #
  # plot(l$rotation[,1,drop=F],p$eigengenes[,"MEblue"])

  p <- WGCNA::moduleEigengenes(expr = expr[, clusters$gene],
                               colors = clusters$module,
                               align = "along average",
                               scale = TRUE)


  p$eigengenes
  p$averageExpr %>% dim

  pp <- firstPC(expr = expr[, clusters$gene],
               colors = clusters$module,
               align = "along average",
               scale = TRUE)

  pp %>% names
  pp$varExplained
  pp$averageExpr
  pp$eigengenes
  pp$PC
  fviz_eig(pp$prcompObj[[6]])

  plot(pp$PC[,2], prcomp(expr[, which(clusters$module %in% "blue")], scale. = TRUE,
       center = TRUE)$x[,2])
  pp$prcompObj[[1]]$rotation[,1]

pp$PC
  pp$averageExpr
  barplot(pp$PC)

  # this plots the eigenvector against the average expression
  # to show the effect of the "along average" argument
  cbind(pp$PC,pp$averageExpr) %>%
    mutate(id = 1:n) %>%
    gather(type, value, -id) %>%
    separate(type, c("type","module")) %>%
    spread(type,value) %>%
    magrittr::set_colnames(c("id","module","average", "PC")) %>%
    ggplot(.,aes(x = average, y = PC)) + geom_point() + facet_grid(~module) +
    theme_bw()

  par(mfrow = c(3,3))
  lapply(pp$prcompObj, get_pca_var)


    p$svdobj[[1]]

  l <- prcomp(t(expr[, which(clusters$module %in% "blue")]), scale. = TRUE,
              center = TRUE)


  l$x %>% dim

  dim(expr)

  ll <- prcomp(expr[, which(clusters$module %in% "blue")], scale. = TRUE,
               center = TRUE)

  expr[, which(clusters$module %in% "blue")] %>% dim
  l$rotation %>% dim

  all.equal(l$sdev, ll$sdev)
  l$sdev %>% sum
  ll$sdev %>% sum
  svd
  ll$rotation %>%  dim
  ll$x %>% dim

  plot(l$rotation[,1],ll$rotation[,1])

  plot(l$rotation[,1,drop=T],p$eigengenes[,"MEblue"])
  all.equal(l$rotation[,1,drop=T],p$eigengenes[,"MEblue"])

  stats:::prcomp.default
  stats:::princomp.default

  prcomp(expr[, which(clusters$module %in% "blue")])$x[,1:2]

  plot(l$rotation[,1,drop=F],p$eigengenes[,"MEblue"])

  plot(p$averageExpr[,"AEblue"], p$eigengenes[,"MEblue"])

  cbind(p$eigengenes,p$averageExpr) %>%
    reshape2::melt() %>%
    # dplyr::select(variable)
    # gsub("AE")
    separate(variable, c("type","module"), sep = "E") %>%
    gather()




  dynamicTreeCut::printFlush

  return(clusters)

}


clusterSimilarity(x = corrX,
                  cutMethod = "dynamic",
                  B = 20) %>% print(nrows = Inf)


clusterSimilarity(x = diffCorr, distanceMethod = "euclidean",
                  cutMethod = "fixed",
                  nClusters = 3,
                  clustMethod = "protoclust") %>% print(nrows = Inf)

protocut()
dist(corrX) %>% as.dist()
as.dist(1 - corrX) %>% str



## ---- data ----
n0 = 60
n1 = 90
n = 150
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
dimnames(TOMX)[[1]] <- dimnames(corrX)[[1]]
dimnames(TOMX)[[2]] <- dimnames(corrX)[[2]]

TOMX0 <- TOMsimilarityFromExpr(X0)
class(TOMX0) <- c("similarity",class(TOMX0))
dimnames(TOMX0)[[1]] <- dimnames(corrX)[[1]]
dimnames(TOMX0)[[2]] <- dimnames(corrX)[[2]]

TOMX1 <- TOMsimilarityFromExpr(X1)
class(TOMX1) <- c("similarity",class(TOMX1))
dimnames(TOMX1)[[1]] <- dimnames(corrX)[[1]]
dimnames(TOMX1)[[2]] <- dimnames(corrX)[[2]]


diffTOM <- abs(TOMX1 - TOMX0)
class(diffTOM) <- c("similarity",class(diffTOM))
dimnames(diffTOM)[[1]] <- dimnames(corrX)[[1]]
dimnames(diffTOM)[[2]] <- dimnames(corrX)[[2]]


alpha <- 1.5
Scorr <- abs(corrX0 + corrX1 - alpha * corrX)
class(Scorr) <- c("similarity", class(Scorr))

Stom <- abs(TOMX0 + TOMX1 - alpha * TOMX)
class(Stom) <- c("similarity", class(Stom))

fisherScore <- fisherZ(n0 = n0,corrX0, n1 = n1, corrX1)


## ---- heat-corr-all ----
hc <- hclust(as.dist(1 - corrX), method = "complete")
plot(corrX, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-corr-e0 ----
hc <- hclust(as.dist(1 - corrX0), method = "complete")
plot(corrX0, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-corr-e1 ----
hc <- hclust(as.dist(1 - corrX1), method = "complete")
plot(corrX1, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

# to reproduce the above you have use this code:
# plot(corrX1, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
#      clustering_distance_rows = as.dist(1 - corrX1),
#      clustering_distance_cols = as.dist(1 - corrX1))


## ---- heat-corr-diff ----
plot(diffCorr, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "complete",
     clustering_distance_rows = dist(diffCorr),
     clustering_distance_cols = dist(diffCorr))

## ---- heat-corr-diff-clust ----
hclustComp <- function(x,k) list(cluster = {
  as.numeric(cutree(hclust(dist(x), method = "complete"), k = k))
})

gapScorr <- cluster::clusGap(diffCorr, FUNcluster = hclustComp, K.max = 10, B = 50)
gapScorr
plot(gapScorr, main = "clusGap(., FUN = protoclust, n.start=20, B= 60)")

lapply(c("firstSEmax", "Tibs2001SEmax", "globalSEmax",
         "firstmax", "globalmax"), function(i) maxSE(f = gapScorr$Tab[, "gap"],
                                                     SE.f = gapScorr$Tab[, "SE.sim"],
                                                     method = i,
                                                     SE.factor = 1))

hc <- hclust(as.dist(diffCorr), method = "complete")

plot(diffCorr, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

plot(diffCorr, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "complete",
     clustering_distance_rows = dist(diffCorr),
     clustering_distance_cols = dist(diffCorr),
     cutree_cols = maxSE(f = gapScorr$Tab[, "gap"],
                         SE.f = gapScorr$Tab[, "SE.sim"],
                         method = "Tibs2001SEmax",
                         SE.factor = 1))




## ---- heat-tom-all ----
hc <- hclust(as.dist(1 - TOMX), method = "complete")
plot(TOMX, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-tom-e0 ----
hc <- hclust(as.dist(1 - TOMX0), method = "complete")
plot(TOMX0, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-tom-e1 ----
hc <- hclust(as.dist(1 - TOMX1), method = "complete")
plot(TOMX1, truemodule = truemodule1, cluster_rows = hc, cluster_cols = hc)

## ---- heat-tom-diff ----
plot(diffTOM, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "complete",
     clustering_distance_rows = dist(diffTOM),
     clustering_distance_cols = dist(diffTOM))

## ---- fishers-zstat ----
plot(fisherScore, truemodule = truemodule1, cluster_rows = T, cluster_cols = T)

## ---- cor-scor ----
plot(Scorr, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "complete")

## ---- cor-scor-clust ----
hclustComp <- function(x,k) list(cluster = {
  as.numeric(cutree(hclust(dist(x), method = "complete"), k = k))
})

gapScorr <- cluster::clusGap(diffCorr, FUNcluster = hclustComp, K.max = 10, B = 50)
gapScorr
plot(gapScorr, main = "clusGap(., FUN = protoclust, n.start=20, B= 60)")

lapply(c("firstSEmax", "Tibs2001SEmax", "globalSEmax",
         "firstmax", "globalmax"), function(i) maxSE(f = gapScorr$Tab[, "gap"],
                                                     SE.f = gapScorr$Tab[, "SE.sim"],
                                                     method = i,
                                                     SE.factor = 1))

hc <- hclust(as.dist(Scorr), method = "complete")



## ---- cor-scor-tom ----
plot(Stom, truemodule = truemodule1, cluster_rows = T, cluster_cols = T,
     clustering_method = "complete")

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
