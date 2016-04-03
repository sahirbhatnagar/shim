##################################
# R source code file for functions used in simulation study
# Git: this is on the eclust repo, simulation branch
# Created by Sahir,  April 2, 2016
# Updated:
# Notes:
# This is based on code from Network analysis book by Horvath
##################################


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
#' @description create a function that takes as input, the number of genes,
#' the true beta vector, the gene expression matrix created
#' from the generate_blocks function
#' and returns a list of data matrix, as well as correlation matrices, TOM matrices,
#' cluster information, training and test data
#' @note this function calls the \code{sim_data} to generate phenotype as a
#' function of the gene expression data. This function also returns other information
#' derived from the simulated data including the test and training sets,
#' the correlation and TOM matrices and the clusters.
#' @return list of (in the following order)
#' \describe{
#' \item{beta_truth}{}
#' \item{distance}{}
#' \item{DT}{data.table of simulated data from the \code{sim_data} function}
#' \item{Y}{}
#' \item{X0}{}
#' \item{X1}{}
#' \item{X_train}{}
#' \item{X_test}{}
#' \item{Y_train}{}
#' \item{Y_test}{}
#' \item{DT_train}{}
#' \item{DT_test}{}
#' \item{S0}{}
#' \item{n_clusters}{}
#' \item{clustered_genes_train}{}
#' \item{clustered_genes_test}{}
#' \item{clusters}{}
#' \item{tom_train_all}{}
#' \item{tom_train_diff}{}
#' \item{tom_train_e1}{}
#' \item{tom_train_e0}{}
#' \item{corr_train_all}{}
#' \item{corr_train_diff}{}
#' \item{corr_train_e1}{}
#' \item{corr_train_e0}{}
#' \item{mse_null}{}
#' }
#'
#' @param p number of genes in design matrix
#' @param X gene expression matrix of size n x p using the \code{generate_blocks}
#' function
#' @param beta true beta coefficient vector
#' @param n total number of subjects
#' @param n0 total number of subjects with E=0
#' @param signal_to_noise_ratio signal to noise ratio, default is 4
#' @param cluster_distance character representing which matrix from the training
#' set that you want to use to cluster the genes. Must be one of the following
#' \itemize{
#' \item corr, corr0, corr1, tom, tom0, tom1, diffcorr, difftom
#' }
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

generate_data <- function(p, X, beta, cluster_distance, n,
                          n0, include_interaction = F,
                          signal_to_noise_ratio = 1) {

  names(beta) <- if (include_interaction) {
    c(paste0("Gene",1:p),"E", paste0("Gene",1:p,":E"))
  } else paste0("Gene",1:p)

  # total true beta vector: this includes all the betas for the genes, then the
  # environment beta, then their interactions if interaction is true.
  # This is used to calculate the model error. This is the same as beta,
  # but in matrix form
  beta_truth <- beta %>% as.matrix()

  # Gene names belonging to the active set
  S0 <- names(beta)[which(beta != 0)]

  n1 <- n - n0

  DT <- sim_data(n = n, n0 = n0, p = p, genes = X,
                 include_interaction = include_interaction,
                 E = c(rep(0,n0),rep(1, n1)),
                 beta = beta,
                 signal_to_noise_ratio = signal_to_noise_ratio) %>%
    as.data.frame()

  Y <- DT[,"Y"] %>% as.matrix

  #remove response from X0 and X1
  X0 <- DT[which(DT$E == 0),-1] %>% as.matrix
  X1 <- DT[which(DT$E == 1),-1] %>% as.matrix

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

  # heatmap data
  genes_e0 <- DT_train[which(DT_train$E == 0),paste0("Gene",1:p)] %>% as.matrix
  genes_e1 <- DT_train[which(DT_train$E == 1),paste0("Gene",1:p)] %>% as.matrix

  tom_train_e0 <- WGCNA::TOMsimilarityFromExpr(genes_e0, corType = "pearson",
                                               power = 6, networkType = "signed",
                                               TOMType = "signed")
  tom_train_e1 <- WGCNA::TOMsimilarityFromExpr(genes_e1, corType = "pearson",
                                               power = 6, networkType = "signed",
                                               TOMType = "signed")
  tom_train_diff <- abs(tom_train_e0 - tom_train_e1)
  tom_train_all <- WGCNA::TOMsimilarityFromExpr(rbind(genes_e0,genes_e1),
                                                corType = "pearson", power = 6,
                                                networkType = "signed",
                                                TOMType = "signed")

  corr_train_e0 <- WGCNA::cor(genes_e0)
  corr_train_e1 <- WGCNA::cor(genes_e1)
  corr_train_diff <- abs(corr_train_e1 - corr_train_e0)
  corr_train_all <- WGCNA::cor(rbind(genes_e0,genes_e1))

  # corScor and Fisher Score matrices
  alpha <- 1.5
  Scorr <- abs(corr_train_e0 + corr_train_e1 - alpha * corr_train_all)
  class(Scorr) <- c("similarity", class(Scorr))

  Stom <- abs(tom_train_e1 + tom_train_e0 - alpha * tom_train_all)
  class(Stom) <- c("similarity", class(Stom))

  fisherScore <- fisherZ(n0 = n0, cor0 = corr_train_e0,
                         n1 = n1, cor1 = corr_train_e1)

  class(tom_train_all) <- append(class(tom_train_all), "similarity")
  class(tom_train_diff) <- append(class(tom_train_diff), "similarity")
  class(tom_train_e1) <- append(class(tom_train_e1), "similarity")
  class(tom_train_e0) <- append(class(tom_train_e0), "similarity")
  class(corr_train_all) <- append(class(corr_train_all), "similarity")
  class(corr_train_diff) <- append(class(corr_train_diff), "similarity")
  class(corr_train_e1) <- append(class(corr_train_e1), "similarity")
  class(corr_train_e0) <- append(class(corr_train_e0), "similarity")

  # Folds for Cross validation
  folds_train <- caret::createFolds(Y_train, k = 10, list = T)
  DT_train_folds <- lapply(folds_train, function(i) DT_train[-i,])
  X_train_folds <- lapply(DT_train_folds, function(i) i[,-grep("Y",colnames(i))])
  Y_train_folds <- lapply(DT_train_folds, function(i) i[,grep("Y",colnames(i))])


  # clusters
  distance <- switch(cluster_distance,
                     corr = corr_train_all,
                     corr0 = corr_train_e0,
                     corr1 = corr_train_e1,
                     diffcorr = corr_train_diff,
                     difftom = tom_train_diff,
                     tom0 = tom_train_e0,
                     tom1 = tom_train_e1,
                     tom = tom_train_all)

  print(paste("The dimension of the distance matrix is",
              dim(distance)[1], "by", dim(distance)[2],
              ". Clutering is done on ",cluster_distance,
              " matrix"))

  print(" starting hierarchical clustering ")

  cl <- hclust(
    as.dist(
      if (cluster_distance %in% c("diffcorr","difftom")) {
        distance
      } else 1 - distance
    ), method = "average"
  )

  cuttree <- dynamicTreeCut::cutreeDynamic(
    cl,
    distM = if (cluster_distance %in% c("diffcorr","difftom")) {
      distance
    } else 1 - distance,
    cutHeight = 0.995, deepSplit = 1,
    pamRespectsDendro = FALSE,
    minClusterSize = 30)

  # check if all cluster groups are 0 which means no cluster
  # assignment and everyone is in their own group
  clusters <- data.table(gene = paste0("Gene",1:p),
                         cluster = if (all(cuttree == 0)) 1:p else cuttree)
  setkey(clusters,"cluster")
  clusters.dat <- table(cuttree)  %>% as.data.frame

  # the cluster assignment of '0' means doesnt belong to a cluster
  n_clusters <- if (all(cuttree == 0)) {
    p
  } else clusters[cluster != 0]$cluster %>% unique %>% length

#   # not being used anywhere
#   clustered_genes_train <- lapply(1:n_clusters, function(i) {
#     DT_train[,clusters[cluster == i]$gene] %>%
#       scale(center = T, scale = F) %>%
#       as.matrix})
#
#   # not being used anywhere
#   clustered_genes_test <- lapply(1:n_clusters, function(i) {
#     DT_test[,clusters[cluster == i]$gene] %>%
#       scale(center = T, scale = F) %>%
#       as.matrix})

  print(paste("Clustering done. The number of clusters is", n_clusters))

  if (include_interaction) {
    gene_groups = copy(clusters)
    gene_groups[,gene := paste0(gene,":E")]
    gene_groups <- rbind(clusters,gene_groups) %>% setkey(cluster)

    pf_temp <- gene_groups[,.N, by = cluster][,pf := sqrt(N)] %>% setkey(cluster)
    gene_groups_inter <- rbind(pf_temp[gene_groups],
                               data.table(cluster = n_clusters + 1, N = 1, pf = 0, gene = "E"))
  }

  # K-Means Clustering on difference in correlation matrices
  # nstart 100 times and then takes the permutation that minimisez the criterion
  #
  #   hier_clust <- hclust(as.dist(tom_train_e1/tom_train_e0), method = "ward.D2")
  #   cutree(hier_clust,2) %>% table
  #   plot(hier_clust)
  #
  #   plot(tom_train_diff)
  #   plot(corr_train_diff)

  #   kmeans_clust <- kmeans(corr_train_diff,2, nstart = 10)
  #   kmeans_clusters <- data.table(gene = paste0("Gene",1:p),cluster = kmeans_clust$cluster)
  #   setkey(kmeans_clusters,"cluster")

  clEclust <- hclust(as.dist(tom_train_diff), method = "average")
  cuttreeEclust <- dynamicTreeCut::cutreeDynamic(clEclust,
                                                 distM = tom_train_diff,
                                                 cutHeight = 0.995, deepSplit = 1,
                                                 pamRespectsDendro = FALSE,
                                                 minClusterSize = 30)
  kmeans_clusters <- data.table(gene = paste0("Gene",1:p), cluster = if (all(cuttreeEclust == 0)) 1:p else cuttreeEclust)
  setkey(kmeans_clusters,"cluster")

  DT <- DT %>% as.matrix
  class(DT) <- append(class(DT),"eset")
  result <- list(beta_truth = beta_truth, distance = distance, DT = DT,
                 Y = Y, X0 = X0, X1 = X1, X_train = X_train, X_test = X_test,
                 Y_train = Y_train, Y_test = Y_test, DT_train = DT_train,
                 DT_test = DT_test, S0 = S0,
                 n_clusters = n_clusters, clustered_genes_train = clustered_genes_train,
                 clustered_genes_test = clustered_genes_test, kmeans_clusters = kmeans_clusters,
                 clusters = clusters, gene_groups_inter = if (include_interaction) gene_groups_inter else NULL,
                 tom_train_all = tom_train_all, tom_train_diff = tom_train_diff, tom_train_e1 = tom_train_e1,tom_train_e0 = tom_train_e0,
                 corr_train_all = corr_train_all, corr_train_diff = corr_train_diff, corr_train_e1 = corr_train_e1, corr_train_e0 = corr_train_e0,
                 mse_null = mse_null, DT_train_folds = DT_train_folds, X_train_folds = X_train_folds, Y_train_folds = Y_train_folds)
  return(result)
}
