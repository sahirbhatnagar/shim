##################################
# R source code file for functions used in simulation study
# Git: this is on the eclust repo, simulation branch
# Created by Sahir,  April 2, 2016
# Updated:
# Notes:
# This is based on code from Network analysis book by Horvath
##################################


## ---- functions ----

"%ni%" <- Negate("%in%")

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
#' @rdname clusterSimilarity
firstPC <- function(expr, colors, exprTest, impute = TRUE, nPC = 1,
                    align = "along average", excludeGrey = FALSE,
                    grey = if (is.numeric(colors)) 0 else "grey",
                    subHubs = TRUE, trapErrors = FALSE,
                    returnValidOnly = trapErrors, softPower = 6, scale = TRUE,
                    verbose = 0, indent = 0) {

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

  PCTest <- data.frame(matrix(NA, nrow = dim(exprTest)[[1]],
                              ncol = length(modlevels)))

  # list to store prcomp objects
  prcompObj <- vector("list", length(modlevels))


  # this is the average expression in a module for each subject
  # so this is a n x length(modlevels) matrix
  averExpr <- data.frame(matrix(NA, nrow = dim(expr)[[1]],
                                ncol = length(modlevels)))

  averExprTest <- data.frame(matrix(NA, nrow = dim(exprTest)[[1]],
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
  names(PC) = paste("pc", modlevels, sep = "")
  names(averExpr) = paste("avg", modlevels, sep = "")
  names(PCTest) = paste("pc", modlevels, sep = "")
  names(averExprTest) = paste("avg", modlevels, sep = "")

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
    datModuleTest <- as.matrix(exprTest[, restrict1])

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

    averExpr[,i] <- rowMeans(datModule, na.rm = TRUE)
    averExprTest[,i] <- rowMeans(datModuleTest, na.rm = TRUE)

    varExpl[[i]] <- factoextra::get_eigenvalue(prcompObj[[i]])[1,"variance.percent"]
    # corAve = cor(averExpr[,i], prcompObj[[i]]$rotation[,1],
    #              use = "p")
    # if (!is.finite(corAve)) corAve = 0
    # if (corAve < 0) prcompObj[[i]]$rotation[,1] = -prcompObj[[i]]$rotation[,1]

    PC[, i] <- predict(prcompObj[[i]])[,1]
    PCTest[, i] <- predict(prcompObj[[i]], newdata = datModuleTest)[,1]
    # plot(PC[, i], prcompObj[[i]]$x[,1])
    #means[i] <- prcompObj[[i]]$center
    #sds[i] <- prcompObj[[i]]$scale


  }

  list(eigengenes = eigenVectors, averageExpr = averExpr,
       averageExprTest = averExprTest,
       varExplained = varExpl, validColors = validColors,
       PC = PC, PCTest = PCTest, prcompObj = prcompObj,
       nclusters = length(modlevels))
}



#' Cluster similarity matrix and return cluster membership of each gene
#' as well as average expression per module, 1st PC of each module, and 1st PC
#' of each module from Test dataset
#' @param x similarity matrix. must have non-NULL dimnames i.e., the rows and
#'   columns should be labelled preferable "Gene1, Gene2, ..."
#' @param expr gene expression data (training set). rows are people, columns are
#'   genes
#' @param exprTest gene expression test set
#' @param distanceMethod  one of "euclidean","maximum","manhattan", "canberra",
#'   "binary","minkowski" to be passed to \code{dist} function. If missing, then
#'   this function will take 1-x as the dissimilarity measure. This
#'   functionality is for diffCorr,diffTOM, fisherScore matrices which need to be
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
                              exprTest,
                              distanceMethod,
                              clustMethod = c("hclust", "protoclust"),
                              cutMethod = c("dynamic","gap", "fixed"),
                              nClusters,
                              method = c("complete", "average", "ward.D2",
                                         "single", "ward.D", "mcquitty",
                                         "median", "centroid"),
                              K.max = 10, B = 50) {

  # x = corrX ; expr = X
  # exprTest = X[sample(seq_len(nrow(X)),nrow(X), replace = TRUE ),]
  # dim(X) ; dim(expr) ; dim(exprTest)
  # clustMethod = c("hclust")
  # cutMethod = c("dynamic")
  # nClusters = 6
  # method = c("complete")
  # summary = c("pca")
  # K.max = 10; B = 50
  # distance = as.dist(1 - x)

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
                                  #cutHeight = 0.995,
                                  deepSplit = 1,
                                  pamRespectsDendro = FALSE,
                                  minClusterSize = 20)
                              } else {
                                hcMod <- hc
                                class(hcMod) <- "hclust"
                                dynamicTreeCut::cutreeDynamic(
                                  hcMod,
                                  distM = as.matrix(distance),
                                  #cutHeight = 0.995,
                                  deepSplit = 1,
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
  # plot(clustAssignment)
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

  # this plots the eigenvector against the average expression
  # to show the effect of the "along average" argument
  # cbind(pp$PC,pp$averageExpr) %>%
  #   mutate(id = 1:n) %>%
  #   gather(type, value, -id) %>%
  #   separate(type, c("type","module")) %>%
  #   spread(type,value) %>%
  #   magrittr::set_colnames(c("id","module","average", "PC")) %>%
  #   ggplot(.,aes(x = average, y = PC)) + geom_point() + facet_grid(~module) +
  #   theme_bw()

  pp <- firstPC(expr = expr[, clusters$gene],
                exprTest = exprTest[, clusters$gene],
                colors = clusters$module,
                align = "along average",
                scale = TRUE)

  # clusters
  # pp %>% names
  # pp$PCTest
  #
  # pp$varExplained
  # pp$averageExpr
  # pp$eigengenes
  # pp$PC


  return(list(clusters = clusters, pcInfo = pp))

}


#' @rdname simulated_data
#' @export
sim_data <- function(n , n0 , p , genes,
                     E, signal_to_noise_ratio = 1,
                     include_interaction = F,
                     beta = NULL) {

  # number of subjects with E=1
  #n1 = n - n0

  # not used for now
  # The coefficients are constructed to have alternating signs and to be exponentially
  # decreasing. (Friedman et al 2010, Journal of Statistical Software)
  # beta <- (-1)^{1:p} * exp(-(2*1:p-1)/100)
  # beta <- (-1)^{1:p} * exp(-(2*1:p - 1 )/600)*sin(1:p)

  #   genes = X
  #   signal_to_noise_ratio = 4
  #   n0 = n1 = 100
  #   E = c(rep(0,n0),rep(1, n1))
  #   beta = c(rnorm(1000),0, rep(0,1000));include_interaction = T
  #   beta = c(rnorm(1000));include_interaction = F

  if (include_interaction) {
    DT <- cbind(genes,E) %>% as.data.table()
    alloc.col(DT,2*p + 1) %>% invisible()
    indx <- grep('Gene', colnames(DT))

    for (j in indx){
      set(DT, i = NULL, j = paste0("Gene",j,":E"), value = DT[[j]]*DT[['E']])
    }
  } else {
    DT <- genes
  }

  y.star <- {DT %>% as.matrix()} %*% beta
  error <- rnorm(n)
  k <- sqrt(var(y.star)/(signal_to_noise_ratio*var(error)))

  y <- y.star + k*error

  result <- if (include_interaction) as.matrix(cbind(y,DT)) else as.matrix(cbind(y,DT,E))
  colnames(result)[1] <- "Y"
  class(result) <- append(class(result), "expression")

  return(result)
}


#' Generate data and test and training sets for simulation study
#'
#' @description create a function that takes as input, the number of genes, the
#'   true beta vector, the gene expression matrix created from the
#'   generate_blocks function and returns a list of data matrix, as well as
#'   correlation matrices, TOM matrices, cluster information, training and test
#'   data
#' @note this function calls the \code{sim_data} to generate phenotype as a
#'   function of the gene expression data. This function also returns other
#'   information derived from the simulated data including the test and training
#'   sets, the correlation and TOM matrices and the clusters.
#' @note the PCs and averages need to be calculated in the fitting functions,
#'   because these will change based on the CV fold
#' @return list of (in the following order) \describe{ \item{beta_truth}{}
#'   \item{distance}{} \item{DT}{data.table of simulated data from the
#'   \code{sim_data} function} \item{Y}{} \item{X0}{} \item{X1}{}
#'   \item{X_train}{} \item{X_test}{} \item{Y_train}{} \item{Y_test}{}
#'   \item{DT_train}{} \item{DT_test}{} \item{S0}{} \item{n_clusters}{}
#'   \item{clustered_genes_train}{} \item{clustered_genes_test}{}
#'   \item{clusters}{} \item{tom_train_all}{} \item{tom_train_diff}{}
#'   \item{tom_train_e1}{} \item{tom_train_e0}{} \item{corr_train_all}{}
#'   \item{corr_train_diff}{} \item{corr_train_e1}{} \item{corr_train_e0}{}
#'   \item{mse_null}{} }
#'
#' @param p number of genes in design matrix
#' @param X gene expression matrix of size n x p using the
#'   \code{generate_blocks} function
#' @param beta true beta coefficient vector
#' @param n total number of subjects
#' @param n0 total number of subjects with E=0
#' @param signal_to_noise_ratio signal to noise ratio, default is 4
#' @inheritParams clusterSimilarity
#' @param cluster_distance character representing which matrix from the training
#'   set that you want to use to cluster the genes. Must be one of the following
#'   \itemize{ \item corr, corr0, corr1, tom, tom0, tom1, diffcorr, difftom,
#'   corScor, tomScor, fisherScore }
#' @param EclustDistance character representing which matrix from the training
#'   set that you want to use to cluster the genes based on the environment.
#'   See \code{cluster_distance} for
#'   avaialble options. Should be different from \code{cluster_distance}. For
#'   example, if \code{cluster_distance=corr} and
#'   \code{EclustDistance=fisherScore}. That is, one should be based on
#'   correlations ignoring the environment, and the other should be based on
#'   correlations accounting for the environment. This function will always
#'   return this add on
#'
#' @examples
#' \dontrun{
#' p = 1000
#' n=200;n0=100
#' beta_genes <- c(runif(50,0.9,1.1),
#'                 runif(50, -1.1,-0.9),
#'                 rep(0,900))
#' # gene expression matrix used in sim_data function
#' X <- mapply(generate_blocks,
#'             rho_E0 = c(-0.70, runif(8, 0.01,0.05), 0.70),
#'             rho_E1 = c(0.70, runif(8, 0.01, 0.05), 0.70),
#'             MoreArgs = list(block_size = 100, n = n, n0 = n0), SIMPLIFY = F) %>%
#'   do.call(cbind, . ) %>%
#'   magrittr::set_colnames(paste0("Gene", 1:1000)) %>%
#'   magrittr::set_rownames(paste0("Subject",1:200))
#'
#' cluster_distance <- "corr"
#' generate_data(p = p, n = n, n0 = n0, X = X, beta_genes = beta_genes, cluster_distance = "corr")
#' }
#' @rdname simulated_data
#' @export

generate_data <- function(p, X, beta,
                          cluster_distance = c("corr", "corr0", "corr1", "tom",
                                               "tom0", "tom1", "diffcorr",
                                               "difftom","corScor", "tomScor",
                                               "fisherScore"),
                          n, n0, include_interaction = F,
                          signal_to_noise_ratio = 1,
                          EclustAddDistance = c("fisherScore", "corScor"),
                          clustMethod = c("hclust", "protoclust"),
                          cutMethod = c("dynamic","gap", "fixed"),
                          distanceMethod = c("euclidean","maximum", "manhattan",
                                             "canberra", "binary", "minkowski"),
                          nClusters,
                          method = c("complete", "average", "ward.D2",
                                     "single", "ward.D", "mcquitty",
                                     "median", "centroid"),
                          K.max = 10, B = 10) {

  # p = p; X = X ; beta = beta
  # n = n; n0 = n0
  # cluster_distance = "corr"
  # include_interaction = T
  # signal_to_noise_ratio = 0.5
  # clustMethod = "hclust" ; cutMethod = "dynamic";method="complete";
  # distanceMethod = "euclidean"
  # EclustAdd = TRUE ; EclustAddDistance = "fisherScore"


  method <- match.arg(method)
  cutMethod <- match.arg(cutMethod)
  clustMethod <- match.arg(clustMethod)
  distanceMethod <- match.arg(distanceMethod)
  cluster_distance <- match.arg(cluster_distance)
  EclustAddDistance <- match.arg(EclustAddDistance)


  names(beta) <- if (include_interaction) {
    c(paste0("Gene",1:p),"E", paste0("Gene",1:p,":E"))
  } else paste0("Gene",1:p)

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

  # tom_train_e0 <- WGCNA::TOMsimilarityFromExpr(genes_e0, corType = "pearson",
  #                                              power = 6, networkType = "signed",
  #                                              TOMType = "signed")
  # tom_train_e1 <- WGCNA::TOMsimilarityFromExpr(genes_e1, corType = "pearson",
  #                                              power = 6, networkType = "signed",
  #                                              TOMType = "signed")
  # tom_train_diff <- abs(tom_train_e0 - tom_train_e1)
  # tom_train_all <- WGCNA::TOMsimilarityFromExpr(genes_all,
  #                                               corType = "pearson", power = 6,
  #                                               networkType = "signed",
  #                                               TOMType = "signed")

  corr_train_e0 <- WGCNA::cor(genes_e0)
  corr_train_e1 <- WGCNA::cor(genes_e1)
  corr_train_diff <- abs(corr_train_e1 - corr_train_e0)
  corr_train_all <- WGCNA::cor(genes_all)

  # corScor and Fisher Score matrices
  alpha <- 1.5
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
  res <- if (cluster_distance %in% c("diffcor","difftom",
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


  resEclust <- if (EclustAddDistance %in% c("diffcor","difftom",
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

  # need to calculate penalty factors for group lasso
  # I put all main effects and interactions of a given module in the same group
  # and the size of the penalty factor is sqrt(size of module), where the
  # size of the module includes both main and interaction effects
  # environment does not get penalized
  if (include_interaction) {

    gene_groups = copy(clustersAll)
    gene_groups[, gene := paste0(gene,":E")]
    gene_groups <- rbind(clustersAll,gene_groups) %>% setkey(cluster)

    pf_temp <- gene_groups[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)

    gene_groups_inter <- rbind(pf_temp[gene_groups],
                               data.table(cluster = n_clusters_All + 1, N = 1,
                                          pf = 0, gene = "E", module = "empty"))

    gene_groups_Addon = copy(clustersAddon)
    gene_groups_Addon[, gene := paste0(gene,":E")]
    gene_groups_Addon <- rbind(clustersAddon, gene_groups_Addon) %>% setkey(cluster)

    pf_temp_Addon <- gene_groups_Addon[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)

    gene_groups_inter_Addon <- rbind(pf_temp_Addon[gene_groups_Addon],
                               data.table(cluster = n_clusters_Addon + 1, N = 1,
                                          pf = 0, gene = "E", module = "empty"))
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
                 # tom_train_all = tom_train_all, tom_train_diff = tom_train_diff,
                 # tom_train_e1 = tom_train_e1,tom_train_e0 = tom_train_e0,
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


# filtering on univariate p-value, taking the lowet XX percent and then running
# multiple linear regression on it
uniFit <- function(train,
                   test,
                   s0,
                   percent = 0.05,
                   stability = F,
                   include_E = F,
                   include_interaction = F,
                   filter_var = F,
                   p = 1000,
                   true_beta) {

  # train: training data, response column should be called "Y"
  #        and the covariates should be labelled Gene1, Gene2,...
  # test: test data, response column should be called "Y"
  #        and the covariates should be labelled Gene1, Gene2,...
  # s0: vector of active betas
  # percent: percentage threshold of highest significant genes
  # note that if youre fitting interactions, the percent should be such that
  # the 2*percent*p is still less than the sample size, else the model wont have any degrees of freedom
  # stability: logical indicating if you just need coefficient values
  # for calculating a stability criterion, or if you need all mse, r2, ect.
  # calculations done
  #         train = result[["DT_train"]] ; test = result[["DT_test"]] ; percent = 0.05 ; stability = F;
  #         include_E = F; include_interaction = F; filter_var=F; p = 1000; s0 = result[["S0"]]
  #         true_beta = result[["beta_truth"]]

  if (include_E == F & include_interaction == T) stop("include_E needs to be T if you want to include interactions")

  # this is used to create univariate filter. If interaction is true, then we filter on the interaction term
  # and the model used to calculate the interaction p-value is Y ~ Gene + E + Gene:E
  # else we just filter on univariate p-value (Y~Gene)
  res <- ldply(paste0("Gene",1:p), function(i) {

    #i="Gene500"
    fit <- lm(as.formula(paste("Y ~",i, if (include_E) "+E", if (include_interaction) paste0("+",i,":E"))),data = train)
    #fit %>% summary
    coef.index <- if (include_interaction) paste0(i,":E") else i

    data.frame(pvalue = as.numeric(pt(abs(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
                                      df = fit$df.residual, lower.tail = F)*2),
               test.stat = as.numeric(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
               'mean' = mean(train[,i]),
               'sd' = sd(train[,i]))
  })

  rownames(res) <- if (include_interaction) paste0("Gene",1:p,":E") else paste0("Gene",1:p)

  top.percent <- percent*p

  uni.S.hat <- if (filter_var) rownames(res[order(res$sd, decreasing = T)[1:top.percent],]) else rownames(res[order(res$pvalue)[1:top.percent],])
  # extract gene names if there is interaction, else just use gene names. this is for the prediction step below
  # and fitting the multivariate model.
  gene.names <- if (include_interaction) stringr::str_sub(string = uni.S.hat, end = -3) else uni.S.hat

  # when fitting the multivariate model with interactions, we must also include main effects
  lm.model <- lm(as.formula(paste("Y ~",paste0(uni.S.hat,collapse = "+"),
                                  if (include_interaction) paste("+",
                                                                 paste0(gene.names, collapse = "+"),"+E"))),
                 data = train)
  #lm.model %>% summary

  lm.oracle <- lm(as.formula(paste("Y ~",paste0(s0,collapse = "+"))), data = train)
  #lm.oracle %>% summary()
  coefs <- if (include_interaction)
    data.table::data.table(Gene = c(paste0("Gene",1:p),"E",rownames(res)),
                           coef.est = rep(0,2*length(rownames(res)) + 1)) else
                             data.table::data.table(Gene = rownames(res),
                                                    coef.est = rep(0,length(rownames(res))))

  coefs[Gene %in% c(uni.S.hat, gene.names, if (include_interaction) "E"), coef.est := coef(lm.model)[-1]]
  #coefs[coef.est!=0]
  #  return(coefs)
  #} else {
  #summary(lm.model)

  if (stability) {
    return(coefs)
  } else {
    lm.pred <- predict.lm(lm.model, newdata = test[,c(gene.names, if (include_interaction) "E")])
    lm.pred.oracle <- predict.lm(lm.oracle, newdata = test[,s0])

    y_test <- test[ , "Y"]

    # Mean Squared Error
    uni.mse <- crossprod(lm.pred - y_test)/length(y_test)
    uni.mse.oracle <- crossprod(lm.pred.oracle - y_test)/length(y_test)

    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)

    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    uni.r2 <- (mse_null - uni.mse)/mse_null

    uni.adj.r2 <- 1 - (1 - uni.r2)*(length(y_test) - 1)/(length(y_test) - length(uni.S.hat) - 1)

    # model error
    identical(true_beta %>% rownames(),coefs[["Gene"]])
    uni.model.error <- {(true_beta - coefs[["coef.est"]]) %>% t} %*% WGCNA::cor(test[,coefs$Gene]) %*% (true_beta - coefs[["coef.est"]])

    # crossprod above gives same result as
    # sum((lm.pred - DT.test[,"V1"])^2)/nrow(DT.test)

    #uni.S.hat %in% S0
    # True Positive Rate
    uni.TPR <- length(intersect(c(gene.names,uni.S.hat), s0))/length(s0)

    # True negatives
    trueNegs <- setdiff(colnames(train), s0)

    # these are the terms which the model identified as zero
    modelIdentifyZero <- setdiff(colnames(train),uni.S.hat)

    # how many of the terms identified by the model as zero, were actually zero
    sum(modelIdentifyZero %in% trueNegs)

    # False Positive Rate = FP/(FP + TN)
    uni.FPR <- sum(uni.S.hat %ni% s0)/(sum(uni.S.hat %ni% s0) + sum(modelIdentifyZero %in% trueNegs))


    ls <- list(uni.mse = as.numeric(uni.mse), uni.r2 = as.numeric(uni.r2),
               uni.adj.r2 = as.numeric(uni.adj.r2), uni.S.hat = length(uni.S.hat),
               uni.TPR = uni.TPR, uni.FPR = uni.FPR, uni.relative.mse = uni.mse/uni.mse.oracle, uni.model.error = uni.model.error)
    names(ls) <- c(paste0("uni",ifelse(filter_var,"_na","_na"),ifelse(include_E,"_lm","_lm"),ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0("uni",ifelse(filter_var,"_na","_na"),ifelse(include_E,"_lm","_lm"),ifelse(include_interaction,"_yes","_no"),"_r2"),
                   paste0("uni",ifelse(filter_var,"_na","_na"),ifelse(include_E,"_lm","_lm"),ifelse(include_interaction,"_yes","_no"),"_adjr2"),
                   paste0("uni",ifelse(filter_var,"_na","_na"),ifelse(include_E,"_lm","_lm"),ifelse(include_interaction,"_yes","_no"),"_Shat"),
                   paste0("uni",ifelse(filter_var,"_na","_na"),ifelse(include_E,"_lm","_lm"),ifelse(include_interaction,"_yes","_no"),"_TPR"),
                   paste0("uni",ifelse(filter_var,"_na","_na"),ifelse(include_E,"_lm","_lm"),ifelse(include_interaction,"_yes","_no"),"_FPR"),
                   paste0("uni",ifelse(filter_var,"_na","_na"),ifelse(include_E,"_lm","_lm"),ifelse(include_interaction,"_yes","_no"),"_relmse"),
                   paste0("uni",ifelse(filter_var,"_na","_na"),ifelse(include_E,"_lm","_lm"),ifelse(include_interaction,"_yes","_no"),"_modelerror"))
    return(ls)
  }
}

clust_fun <- function(x_train,
                      x_test,
                      y_train,
                      y_test,
                      s0,
                      summary,
                      model,
                      gene_groups,
                      true_beta = NULL,
                      topgenes = NULL,
                      stability = F,
                      filter = F,
                      include_E = F,
                      include_interaction = F,
                      p = 1000,
                      filter_var = F,
                      clust_type = c("clust","Eclust","Addon")){

# result[["clustersAddon"]] %>% print(nrows=Inf)
#     result %>% names
#                 stability = F; gene_groups = result[["clustersAll"]];
#                 x_train = result[["X_train"]] ; x_test = result[["X_test"]] ;
#                 y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
#                 filter = F; filter_var = F; include_E = T; include_interaction = T;
#                 s0 = result[["S0"]]; p = p ; k_means = NULL
#                 model = "shim"; summary = "pc"; topgenes = NULL;
  #
  #             stability = F; gene_groups = result[["clusters"]][gene %in% result[["S0"]]];
  #             x_train = result[["X_train"]][,result[["S0"]]] ; x_test = result[["X_test"]][,result[["S0"]]] ;
  #             y_train = result[["Y_train"]] ; y_test = result[["Y_test"]];
  #             filter = F; filter_var = F; include_E = F; include_interaction = F;
  #             s0 = result[["S0"]]; p = 1000 ;
  #             model = "lm"; summary = "spc"; topgenes = NULL;

  #true_beta = result[["beta_truth"]]

  # model: lm, lasso, scad, mcp, elasticnet -
  # summary: summary measure either "pc", "avg", "spc"
  # gene_groups: data table with columns 'gene' and 'cluster'
  # filter: T or F based on univariate filter

  clust_type <- match.arg(clust_type)

  message(paste(summary,model,"filter = ", filter, "filter_var = ", filter_var,
              "include_E = ", include_E, "include_interaction = ",
              include_interaction, sep = ","))

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

  # train data which includes the relevant (filtered or not filtered genes and E or not E)
  x_train_mod <- if (filter & !include_E) {x_train[, topgenes] %>% as.data.frame} else if (!filter & include_E) {
    x_train %>% as.data.frame } else if (!filter & !include_E) {
      x_train[,which(colnames(x_train) %ni% "E")] %>% as.data.frame} else if (filter & include_E) {
        x_train[, c(topgenes,"E")] %>% as.data.frame
      }

  # test data
  x_test_mod = if (filter & !include_E) {x_test[, topgenes] %>% as.data.frame} else if (!filter & include_E) {
    x_test %>% as.data.frame } else if (!filter & !include_E) {
      x_test[,which(colnames(x_test) %ni% "E")] %>% as.data.frame} else if (filter & include_E) {
        x_test[, c(topgenes,"E")] %>% as.data.frame
      }

  # this includes interaction terms
  gene.names <- colnames(x_train_mod)[which(colnames(x_train_mod) %ni% "E")]

  # these are only derived on the main effects genes.. E is only included in the model
  PC_and_avg <- firstPC(expr = x_train_mod[,gene_groups$gene],
                        colors = gene_groups$cluster,
                        exprTest = x_test_mod[,gene_groups$gene])

  n.clusters <- PC_and_avg$nclusters

  # the number of clusters in the modified set
  # remove genes which have a cluster assignment of 0, meaning they dont cluster well
  # n.clusters <- PC_and_avg$nclusters
  # cluster_names <- PC_and_avg$validColors
  #
  #
  # clustered.genes <- lapply(cluster_names, function(i) {
  #   x_train_mod[,intersect(gene_groups[cluster == i]$gene, gene.names), drop = FALSE] %>%
  #     scale(center = T, scale = F) %>%
  #     as.matrix
  # })
  #
  # # remove clusters that have no data in them due to filtering
  # clustered.genes <- Filter(function(i) ncol(i) > 0, clustered.genes)
  # print(sapply(clustered.genes, dim))

  # this contains either the averages or PCs for each module in a data.frame
  clust_data <- switch(summary,
                       avg = PC_and_avg$averageExpr,
                       pc = PC_and_avg$PC,
                       spc = {
                         # works but never used
                         out <- lapply(clustered.genes,
                                       function(i) PMA::SPC.cv(i,
                                                               sumabsvs =  seq(1, sqrt(ncol(i)), len = 10),
                                                               nfolds = 5, niter = 5, center = F)

                         )

                         out_best_lambda <- lapply(out, function(i) i$bestsumabsv)

                         out_best <- mapply(PMA::SPC,
                                            x = clustered.genes,
                                            sumabsv = out_best_lambda,
                                            MoreArgs = list(K = 1, center = F),
                                            SIMPLIFY = F)

                         out_best_v <- mapply(function(spc, names) {
                           spc$v %>% magrittr::set_rownames(colnames(names))
                         }, spc = out_best, names = clustered.genes,SIMPLIFY = F)

                         # output the sparse factors as well as model fit for use with testing below
                         list(mapply(new.pc, data = clustered.genes,
                                     loadings = out_best_v), out_best, out_best_v)
                       })

  # if (summary == "spc") {
  #   out_best <- clust_data[[2]]
  #   out_best_v <- clust_data[[3]] # used to calculate S.hat (number of non zero components)
  #
  #   # these are the gene names corresponding to the non zero loadings from sparsePCA
  #   non_zero_pc_components <- lapply(out_best_v, function(i)
  #     i[which(i != 0), ,drop = F] %>% rownames) %>% unlist
  # }
  #
  # clust_data <- if (summary == "spc") clust_data[[1]] else clust_data

  # colnames(clust_data) <- paste0(summary,cluster_names)

  # head(clust_data)

  ml.formula <- if (include_interaction & include_E) {
    paste0("y_train ~","(",paste0(colnames(clust_data), collapse = "+"),")*E") %>% as.formula} else if (!include_interaction & include_E) {
      paste0("y_train ~",paste0(colnames(clust_data), collapse = "+"),"+E") %>% as.formula} else if (!include_interaction & !include_E) {
        paste0("y_train ~",paste0(colnames(clust_data), collapse = "+")) %>% as.formula
      }

  # this is the same as ml.formula, except without the response.. this is used for
  # functions that have the x = and y = input instead of a formula input
  model.formula <- if (include_interaction & include_E) {
    paste0("~ 0+(",paste0(colnames(clust_data), collapse = "+"),")*E") %>% as.formula} else if (!include_interaction & include_E) {
      paste0("~0+",paste0(colnames(clust_data), collapse = "+"),"+E") %>% as.formula} else if (!include_interaction & !include_E) {
        paste0("~0+",paste0(colnames(clust_data), collapse = "+")) %>% as.formula
      }

  X.model.formula <- model.matrix(model.formula, data = if (include_E) cbind(clust_data,x_train_mod[,"E", drop = F]) else
    clust_data %>% as.data.frame )
  df <- X.model.formula %>% cbind(y_train) %>% as.data.frame()

  clust_train_model <- switch(model,
                              lm = lm(ml.formula , data = df),
                              lasso = {if (n.clusters != 1) {
                                cv.glmnet(x = X.model.formula, y = y_train, alpha = 1)
                              } else NA },
                              elasticnet = {if (n.clusters != 1) {
                                cv.glmnet(x = X.model.formula, y = y_train, alpha = 0.5)
                              } else NA },
                              scad = {
                                cv.ncvreg(X = X.model.formula, y = y_train,
                                          family = "gaussian", penalty = "SCAD")
                              },
                              mcp = {
                                cv.ncvreg(X = X.model.formula, y = y_train,
                                          family = "gaussian", penalty = "MCP")
                              },
                              shim = {
                                require(doMC)
                                registerDoMC(cores = 2)
                                cv.shim(x = X.model.formula, y = y_train,
                                        main.effect.names = c(colnames(clust_data), if (include_E) "E"),
                                        interaction.names = setdiff(colnames(X.model.formula),c(colnames(clust_data),"E")),
                                        max.iter = 500, initialization.type = "univariate",
                                        verbose = FALSE, threshold = 1e-06,
                                        nlambda.gamma = 10, nlambda.beta = 10, nlambda = 100)
                              })

  # here we give the coefficient stability on the clusters and not the individual genes
  coefs <- switch(model,
                  lm = data.table::data.table(Gene = names(clust_train_model$coefficients),
                                              coef.est = coef(clust_train_model)),
                  lasso = {
                    # need to return all 0's if there is only 1 cluster since lasso wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {
                      coef(clust_train_model, s = "lambda.min") %>%
                        as.matrix %>%
                        as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  },
                  elasticnet = {
                    # need to return all 0's if there is only 1 cluster since lasso wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {
                      coef(clust_train_model, s = "lambda.min") %>%
                        as.matrix %>%
                        as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  },
                  scad = {
                    # need to return all 0's if there is only 1 cluster since lasso wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {

                      coef(clust_train_model, lambda = clust_train_model$lambda.min) %>%
                        as.matrix %>%
                        as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  },
                  mcp = {
                    # need to return all 0's if there is only 1 cluster since lasso wont run with only 1 predictor
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {

                      coef(clust_train_model, lambda = clust_train_model$lambda.min) %>%
                        as.matrix %>%
                        as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  },
                  shim = {
                    dat <- data.table::data.table(Gene = colnames(X.model.formula),
                                                  coef.est = rep(0, ncol(X.model.formula)))
                    if (n.clusters != 1) {

                      coef(clust_train_model, s = "lambda.min") %>%
                        as.matrix %>%
                        as.data.table(keep.rownames = TRUE) %>%
                        magrittr::set_colnames(c("Gene","coef.est"))
                    } else dat
                  }
  )

  coefs

  if (stability) {
    # remove intercept for stability measures
    return(coefs %>% magrittr::extract(-1, , drop = F))
  } else {

    non_zero_clusters <- coefs[-1, , ][coef.est != 0] %>%
      magrittr::use_series("Gene")

    # non_zero_clusters <- coefs[-1, , ] %>%
    #   magrittr::use_series("Gene")
    #
    # non_zero_clusters <- non_zero_clusters[c(-8,-9)]

    # need to determine which of non_zero_cluters are main effects and which
    # are interactions
    non_zero_clusters_interactions <- grep(":",non_zero_clusters, value = T)
    non_zero_environment <- grep("^E", non_zero_clusters, value = T, ignore.case = TRUE)
    non_zero_clusters_main_effects <- setdiff(non_zero_clusters, c(non_zero_clusters_interactions,
                                                                   non_zero_environment))
    n.non_zero_clusters <- coefs[-1, , ][coef.est != 0] %>%
      magrittr::use_series("Gene") %>%
      length

    # need to get the genes corresponding to the non-zero clusters
    # NOTE: this also includes non-zero cluster:Environment interactions

    # genes corresponding to non-zero clusters
    # return non zero sparse loadings if spc combined with lm, else
    # return non-zero clusters for penalization methods, even if your using spc
    # clust.S.hat <- if (summary == "spc" & model == "lm") non_zero_pc_components else {
    #   gene_groups[cluster %in%
    #                 as.numeric(unlist(stringr::str_extract_all(non_zero_clusters, "\\d+"))),gene] }

    clust.S.hat.main <- gene_groups[cluster %in% as.numeric(unlist(stringr::str_extract_all(non_zero_clusters, "\\d+"))),gene]

    # this is the same as gene_groups, but the gene names contain E
    # so that we can extract the interactions corresponding to the chose clusters
    gene_groups_environment <- copy(gene_groups)
    gene_groups_environment[,gene:=paste0(gene,":E")]

    clust.S.hat.interaction <- gene_groups_environment[cluster %in% as.numeric(unlist(stringr::str_extract_all(non_zero_clusters_interactions, "\\d+"))),gene]

    # this represents all the genes corresponding to the non-zero PC or avg
    clust.S.hat <- c(clust.S.hat.main, non_zero_environment, clust.S.hat.interaction)
    # clustered.genes_test <- lapply(cluster_names, function(i) {
    #   x_test_mod[,intersect(gene_groups[cluster == i]$gene, gene.names), drop = FALSE] %>% scale(center = T, scale = F) %>% as.matrix
    # }
    # )
    #
    # # remove clusters that have no data in them due to filtering
    # clustered.genes_test <- Filter(function(i) ncol(i) > 0, clustered.genes_test)

    clust_data_test <- switch(summary,
                              avg = PC_and_avg$averageExprTest,
                              pc = PC_and_avg$PCTest,
                              spc = {
                                V <- lapply(out_best, function(i) i$v)
                                mapply(new.pc, data = clustered.genes_test,
                                       loadings = V)
                              }
    )

    #colnames(clust_data_test) <- paste0(summary,cluster_names)

    # need intercept for prediction
    model.formula_test <- if (include_interaction & include_E) {
      paste0("~ 1+(",paste0(colnames(clust_data_test), collapse = "+"),")*E") %>% as.formula} else if (!include_interaction & include_E) {
        paste0("~1+",paste0(colnames(clust_data_test), collapse = "+"),"+E") %>% as.formula} else if (!include_interaction & !include_E) {
          paste0("~1+",paste0(colnames(clust_data_test), collapse = "+")) %>% as.formula
        }

    X.model.formula_test <- model.matrix(model.formula_test,
                                         data = if (include_E) {
                                           cbind(clust_data_test,x_test_mod[,"E", drop = F])
                                           } else clust_data_test %>% as.data.frame)



    # True Positive Rate
    clust.TPR <- length(intersect(clust.S.hat, s0))/length(s0)

    # True negatives
    trueNegs <- setdiff(colnames(x_train), s0)

    # these are the terms which the model identified as zero
    modelIdentifyZero <- setdiff(colnames(x_train),clust.S.hat)

    # how many of the terms identified by the model as zero, were actually zero
    sum(modelIdentifyZero %in% trueNegs)

    # False Positive Rate = FP/(FP + TN)
    clust.FPR <- sum(clust.S.hat %ni% s0)/(sum(clust.S.hat %ni% s0) + sum(modelIdentifyZero %in% trueNegs))

    # Mean Squared Error
    clust.mse <- crossprod(X.model.formula_test %*% coefs$coef.est - y_test)/length(y_test)

    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)

    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    clust.r2 <- (mse_null - clust.mse)/mse_null

    clust.adj.r2 <- 1 - (1 - clust.r2)*(nrow(x_test) - 1)/(nrow(x_test) - n.non_zero_clusters - 1)


    ls <- list(clust.mse = clust.mse, clust.r2 = clust.r2, clust.adj.r2 = clust.adj.r2, clust.S.hat = length(clust.S.hat),
               clust.TPR = clust.TPR, clust.FPR = clust.FPR)


    names(ls) <- c(paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_r2"),
                   paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_adjr2"),
                   paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                   paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                   paste0(clust_type,"_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"))
    return(ls)

  }
}



