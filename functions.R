##################################
# R source code file for functions used in Co-expression simulation study
# 
# Created by Sahir, July 20th, 2015
# Updated: 
# NOTE: adding noise to correlation structures based on the paper
# "A method for generating realistic correlation matrices" by Hardin et al.
# Annals of applied stats 2013
##################################


## ---- toeplitz ----

simcorTop <- function(k=6,size=c(10,5,8,2,15,50),rho=c(.7,.7,.5,.9,.85,.4), eidim=3,
                      c=1, addNoise = T) {
  
  # this function calculates an AR(1) Toeplitz matrix
  # with a block diagonal structure. The size and base correlation for each
  # block is user specified
  # k: is the number of groups
  # size: is a vector of length k specifying the size of each group 
  # rho: is a vector of length k specifying base correlation values
  # eidim is the space from which the noise is generated, the smaller the more noise
  # addNoise: logical indicating wether to add noise as in Hardin et al. or not
  # c: parameter to make correlations between genes stronger or weaker
  #    only use if youre not adding noise, else the correlation matrix
  #    wont be positive definite
  
  #epsilon <- (1-max(rho))/(1+max(rho))
  epsilon <- 0.05
  ndim <- sum(size) # dim of correlation matrix
  bigcor<- matrix(rep(0, ndim*ndim), ncol=ndim)
  
  for (i in 1:k){
    
    sigma = 1
    times = 1:size[i]
    
    # simulate Toeplitz like structure
    H <- c*abs(outer(times, times, "-"))
    cor <- sigma * rho[i]^H
    
    if (i==1){bigcor[1:size[1], 1:size[1]] <- cor}
    if (i!=1){bigcor[(sum(size[1:(i-1)]) + 1):sum(size[1:i]),
                     (sum(size[1:(i-1)]) + 1):sum(size[1:i])] <- cor}
  }
  
  class(bigcor) <- append(class(bigcor),"truecorr")
  
  
  if(addNoise) {
    
    diag(bigcor) <- 1 - epsilon
    
    #return(corr.nz)
    
    ### adding noise to the correlation matrix
    
    eivect <- c( )
    for (i in 1:ndim) {
      ei <- runif(eidim, -1, 1)
      eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2) ) )
    }
    
    
    bigE <- t(eivect) %*% eivect
    cor.nz <- bigcor + bigE
    class(cor.nz) <- append(class(cor.nz),"truecorr")
  }
  
  if(addNoise) return(cor.nz) else return(bigcor)
  
}




## ---- hub ----
#  Simulating the Hub Matrix (entries filled in using Toeplitz structure)
# this function calculates a Toeplitz matrix with values descending
# from a user specified maximum to minimum.  The matrix has a 
# block diagonal structure.  The size and base correlation for each
# block is user specified.
# k is the number of groups
# size is a vector of length k specifying the size of each group 
# rho is a vector of length k specifying base correlation values (rho max and rho min)
# rho[1] is the max and rho[2] is the min correlation
# epsilon <- (1-min(rho) - 0.75*min(tau) ) - .01
# tau_k = (max(rho_k) - min(rho_k) )/ (size_k -2) 
# eidim is the space from which the noise is generated, the smaller the more noise
# power = 2 makes the correlations stay high
# power = 0.5 makes the correlations descent rapidly

# rho.func is needed for filling in the rest of the structure of the Hub
# correlation matrix
# r.max is the maximum user specified correlation
# r.min is the minimum user specified correlation
# power is the power at which the correlations descend
# p is the size of the correlation block

rho.func <- function(r.max, r.min, power,p){
  rhovec <-c()
  
  rhovec[1] <- 1
  for(i in 2:p){
    rhovec[i] <- r.max - ((i-2)/(p-2))^power*(r.max-r.min)
  }
  rhovec}

simcor.H <- function(k=5, size=rep(200,5), 
                     rho=rbind(c(.9,.7), c(.7,.7), c(.7,.2), c(.5,.3), c(.9,.85)), power=1,
                     eidim=2, epsilon){
  
  tau_k = sapply(1:k, function(i) (rho[i,1] - rho[i,2])/ (size[i] - 2) )
  #epsilon <- min(sapply(1:k, function(i) {1 - rho[i,1] - 0.75*tau_k[i]}))-0.01
  #epsilon <- 0.29
  
  ndim <- sum(size)# dim of correlation matrix
  bigcor<- matrix(rep(0, ndim*ndim), ncol=ndim)
  
  ### generating the basic correlation matrix
  
  
  for (i in 1:(k) ){
    
    cor <- toeplitz(rho.func(rho[i,1],rho[i,2],power,size[i]) )
    
    if (i==1){bigcor[1:size[1], 1:size[1]] <- cor}
    if (i!=1){bigcor[(sum(size[1:(i-1)]) + 1):sum(size[1:i]),
                     (sum(size[1:(i-1)]) + 1):sum(size[1:i])] <- cor}
  }
  diag(bigcor) <- 1 - epsilon
  
  
  ### adding noise to the correlation matrix
  
  eivect <- c( )
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei / sqrt(sum(ei^2) ) )
  }
  
  
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  class(cor.nz) <- append(class(cor.nz),"truecorr")
  return(cor.nz)
}



## ---- constant ----

#  Simulating the Constant Correlation

# this function simulates a block correlation matrix
# The size and base correlation for each block is user specified
# There is an additional delta parameter for the off diagonal correlations

# k is the number of groups
# size is a vector of length k specifying the size of each group 
# rho is a vector of length k specifying base correlation values
# epsilon <- 0.99 - max(rho)
# eidim is the space from which the noise is generated, the smaller the more noise
# delta is the correlation of the off diagonal blocks

simcor.const = function (k = 5, size = rep(200,5), 
                   rho = c(0.7, 0.7, 0.5, 0.9, 0.85), 
                   delta = 0.39, eidim = 2, epsilon=0.01) 
{
  #epsilon = 0.98 - max(rho)
  #epsilon = 0.0
  ndim <- sum(size)
  bigcor <- matrix(rep(delta, ndim * ndim), ncol = ndim)
  for (i in 1:k) {
    cor <- matrix(rep(rho[i], size[i] * size[i]), ncol = size[i])
    if (i == 1) {bigcor[1:size[1], 1:size[1]] <- cor}
    if (i != 1) {bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]), 
                        (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor}
  }
  diag(bigcor) <- 1 - epsilon
  
  eivect <- c()
  for (i in 1:ndim) {
    ei <- runif(eidim, -1, 1)
    eivect <- cbind(eivect, sqrt(epsilon) * ei/sqrt(sum(ei^2)))
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  class(cor.nz) <- append(class(cor.nz),"truecorr")
  return(cor.nz)
  
}




## ---- simulate-data ----
#' Simulate data with differential correlation depending on a binary environment variable
#' 
#' @description 
#' Generate gene expression data for a block of genes for n subjects.
#' The data is generated depending on the binary environment variable E. The data is
#' generated such that the first n0 rows are subjects who have E = 0, and the next
#' n-n0 rows are subjects with E = 1
#' @name simulated_data
NULL
#' 
#' @return simulated gene expression data 
#' 
#' @param block_size the number of genes in this block
#' @param rho_E0 correlation between -1 and 1 of the genes in the block for 
#' subjects with E = 0
#' @param rho_E1 correlation between -1 and 1 of the genes in the block for 
#' subjects with E = 1
#' @param beta true beta coefficient vector of length 2p+1 if \code{include_interaction} is TRUE. 
#' Must be ordered in this way: first the coefficients for Gene1, Gene2, ..., Genep, E, Gene1:E, 
#' Gene2:E, ..., Genep:E. If \code{include_interaction} is FALSE should be vector of length p.
#' @param n total number of subjects
#' @param n0 number of subjects with E = 0
#' @param genes gene expression data matrix where rows are subjects and columns 
#' are genes. Rownames should be "Subject1", "Subject2", ... and colnames should
#' be "Gene1", "Gene2",... This information is used for annotation heatmaps.
#' @param signal_to_noise_ratio signal to noise ratio see ESL book for details
#' @param include_interaction logical indicating if youre doing an analysis with interactions
#' this will affect all downstream analysis
#' @param E numeric vector for environment status of subjects in \code{genes} 
#' data. Should be 0 for E=0 and 1 for E=1.
#' @note to simulate data using you need to first run the \code{generate_blocks} function
#' to create the \code{genes} and then subsequently pass this to the \code{sim_data} function.
#' \code{sim_data} 
#' @rdname simulated_data
#' @export

generate_blocks <- function(block_size, rho_E0, rho_E1, n, n0) {
  #rho_E0 = -0.55; rho_E1 = 0.75 ; n = 200; n0 = 100 ; block_size = 100
  
  t_E0 <- replicate(block_size, rnorm(n0, sd = 1 - abs(rho_E0))) 
  z_E0 <- replicate(1, rnorm(n0, sd = abs(rho_E0)))  
  
  # this returns an n0 x block_size matrix (subjects are rows, columns are genes)
  # of gene expression
  genes0 <- lapply(1:ncol(t_E0), function(i) { 
    if (rho_E0 < 0) {
      if (i <= block_size/2) t_E0[,i] + z_E0 else t_E0[,i] - z_E0   
    } else {
      t_E0[,i] + z_E0
    }
  }
  ) %>% 
    do.call(cbind, .)
  
  #   genes0 %>% dim
  #   
  #   genes0 %>% cor %>%  pheatmap::pheatmap(cluster_rows = F,
  #                                          cluster_cols = F,
  #                                          show_colnames = F,
  #                                          show_rownames = F)
  #   genes0 %>%
  #     magrittr::set_colnames(paste0("Gene", 1:100)) %>% 
  #     magrittr::set_rownames(paste0("Subject",1:100)) %>% 
  #     pheatmap::pheatmap(cluster_rows = F,
  #                       cluster_cols = F,
  #                       show_colnames = T,
  #                       show_rownames = T,
  #                       color = colorRampPalette((brewer.pal(n = 7, name =
  #                                                                 "Reds")))(100))
  
  t_E1 <- replicate(block_size, rnorm(n - n0, sd = 1 - abs(rho_E1))) 
  z_E1 <- replicate(1, rnorm(n - n0, sd = abs(rho_E1)))  
  
  # this returns an n-n0 x block_size matrix (subjects are rows, columns are genes)
  # of gene expression
  genes1 <- lapply(1:ncol(t_E1), function(i) { 
    if (rho_E1 < 0) {
      if (i <= block_size/2) t_E1[,i] + z_E1 else t_E1[,i] - z_E1   
    } else {
      t_E1[,i] + z_E1
    }
  }
  ) %>% 
    do.call(cbind, .)
  
  #   genes1 %>% cor %>%  pheatmap::pheatmap(cluster_rows = F,
  #                                          cluster_cols = F,
  #                                          show_colnames = F,
  #                                          show_rownames = F)
  #   
  #   genes1 %>%
  #     magrittr::set_colnames(paste0("Gene", 1:100)) %>% 
  #     magrittr::set_rownames(paste0("Subject",101:200)) %>% 
  #     pheatmap::pheatmap(cluster_rows = F,
  #                        cluster_cols = F,
  #                        show_colnames = T,
  #                        show_rownames = T,
  #                        color = colorRampPalette(brewer.pal(n = 7, name =
  #                                                               "Reds"))(100))
  
  return(rbind(genes0,genes1))
  
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
  
#   rm(list = ls())
#   #dev.off()
#   source("/home/data1/share/sy/phd/coexpression/functions.R")
#   # gene expression matrix used in sim_data function
#   p = 1000; n = 400; n0 = 200;n1 = n-n0
#   
#   rho.0.hub <- rbind(c(.05,.01), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5))
#   rho.1.hub <- rbind(c(.9,.7), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5), c(.5,.5))
#   
#   V0 = simcor.H(k = 10, size = rep(100,10), 
#                 rho = rho.0.hub, power=1,
#                 eidim = 2, epsilon = 0.05)
#   V1 = simcor.H(k = 10, size = rep(100,10), 
#                 rho = rho.1.hub, power=1,
#                 eidim = 2, epsilon = 0.05)
#   
#   genes0 <- MASS::mvrnorm(n = n0, mu = rep(0,p), Sigma = V0)
#   genes1 <- MASS::mvrnorm(n = n1, mu = rep(0,p), Sigma = V1)
#   
#   X <- rbind(genes0,genes1) %>%  
#     magrittr::set_colnames(paste0("Gene", 1:p)) %>% 
#     magrittr::set_rownames(paste0("Subject",1:n))
#   
#   
#   # dataset without interaction
#   beta <- c(rep(0,50),runif(50, 0.9,1.1), rep(0,900));include_interaction = F
#   cluster_distance = "corr";signal_to_noise_ratio = 1
    
    
  #============================================#=====================================================#
  
  names(beta) <- if (include_interaction) c(paste0("Gene",1:p),"E",paste0("Gene",1:p,":E")) else paste0("Gene",1:p)
  
  # total true beta vector: this includes all the betas for the genes, then the 
  # environment beta, then their interactions if interaction is true. 
  #  This is used to calculate the model error. This is the same as beta, but in matrix form
  beta_truth <- beta %>% as.matrix()
  
  # Gene names belonging to the active set
  S0 <- names(beta)[which(beta != 0)]
  
  n1 <- n - n0
  
  DT <- sim_data(n = n, n0 = n0, p = p, genes = X,
                 include_interaction = include_interaction,
                 E = c(rep(0,n0),rep(1, n1)),
                 beta = beta,
                 signal_to_noise_ratio = signal_to_noise_ratio) %>% as.data.frame()
  
  Y <- DT[,"Y"] %>% as.matrix
  X0 <- DT[which(DT$E == 0),c(-1,-2)] %>% as.matrix
  X1 <- DT[which(DT$E == 1),c(-1,-2)] %>% as.matrix
  
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
  
  tom_train_e0 <- WGCNA::TOMsimilarityFromExpr(genes_e0, corType = "pearson", power = 6, networkType = "signed",  TOMType = "signed")
  tom_train_e1 <- WGCNA::TOMsimilarityFromExpr(genes_e1, corType = "pearson", power = 6, networkType = "signed",  TOMType = "signed")
  tom_train_diff <- abs(tom_train_e0 - tom_train_e1)
  tom_train_all <- WGCNA::TOMsimilarityFromExpr(rbind(genes_e0,genes_e1),  corType = "pearson", power = 6, networkType = "signed",  TOMType = "signed")
  
  corr_train_e0 <- WGCNA::cor(genes_e0)
  corr_train_e1 <- WGCNA::cor(genes_e1)
  corr_train_diff <- abs(corr_train_e0 - corr_train_e1)
  corr_train_all <- WGCNA::cor(rbind(genes_e0,genes_e1))
  
#   for (i in c(tom_train_all,tom_train_diff,tom_train_e1,tom_train_e0,
#               corr_train_all,corr_train_diff,corr_train_e1,corr_train_e0)) {
#     class(i) <- append(class(i), "similarity")
#   }
  
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
  
  print(paste("The dimension of the distance matrix is", dim(distance)[1], "by", dim(distance)[2],
              ". Clutering is done on ",cluster_distance, " matrix"))
  
  print(" starting hierarchical clustering ")
  
  cl <- hclust(as.dist(if (cluster_distance %in% c("diffcorr","difftom")) distance else 1 - distance), method = "ward.D2")
  cuttree <- dynamicTreeCut::cutreeDynamic(cl, 
                                           distM = if (cluster_distance %in% c("diffcorr","difftom")) distance else 1 - distance, 
                                           deepSplit = 1)
  # check if all cluster groups are 0 which means no cluster assignment and everyone is in their own group
  clusters <- data.table(gene = paste0("Gene",1:p),cluster = if (all(cuttree == 0)) 1:p else cuttree)
  setkey(clusters,"cluster")
  clusters.dat <- table(cuttree)  %>% as.data.frame
  
  # the cluster assignment of '0' means doesnt belong to a cluster
  n_clusters <- if (all(cuttree == 0)) p else clusters[cluster != 0]$cluster %>% unique %>% length
  
  clustered_genes_train <- lapply(1:n_clusters, function(i) {
    DT_train[,clusters[cluster == i]$gene] %>% scale(center = T, scale = F) %>% as.matrix})
  
  clustered_genes_test <- lapply(1:n_clusters, function(i) {
    DT_test[,clusters[cluster == i]$gene] %>% scale(center = T, scale = F) %>% as.matrix})
  
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
  
  kmeans_clust <- kmeans(tom_train_diff,2, nstart = 10)
  kmeans_clusters <- data.table(gene = paste0("Gene",1:p),cluster = kmeans_clust$cluster)
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


# this one is not used anymore
# sim.data.fun_revised <- function(n = 200, n0 = 100, p = 1000, 
#                                  signal.to.noise.ratio = 4, 
#                                  beta.genes = rep(1, 1000)){
#   
#   #n = 200; n0 = 100; p = 1000; 
#   #signal.to.noise.ratio = 4;  beta.genes = rep(1, 1000)
#   # n: number of subjects
#   # n0: number of subjects with E=0
#   # p: number of genes to simulate (needs to match the dimension of V0 and V1)
#   # V0: correlation matrix for subjects with E=0
#   # V1: correlation matrix for subjects with E=1
#   # signal.to.noise.ratio: signal to noise ratio see ESL book for details
#   # beta: coefficients for each gene
#   
#   # number of subjects with E=1
#   n1 = n - n0 
#   
#   # covariate vector
#   E = c(rep(0,n0),rep(1, n1))
#   
#   # rows are people, columns are genes
#   genesq0 <- MASS::mvrnorm(n = n0, mu = rep(0,p/4), Sigma = q0)
#   genesq1 <- MASS::mvrnorm(n = n1, mu = rep(0,p/4), Sigma = q1)
#   genesr0 <- MASS::mvrnorm(n = n0, mu = rep(0,p/4), Sigma = r0)
#   genesr1 <- MASS::mvrnorm(n = n1, mu = rep(0,p/4), Sigma = r1)
#   geness0 <- MASS::mvrnorm(n = n0, mu = rep(0,p/4), Sigma = s0)
#   geness1 <- MASS::mvrnorm(n = n1, mu = rep(0,p/4), Sigma = s1)
#   genest0 <- MASS::mvrnorm(n = n0, mu = rep(0,p/4), Sigma = t0)
#   genest1 <- MASS::mvrnorm(n = n1, mu = rep(0,p/4), Sigma = t1)
#   
#   genes <- rbind(cbind(genesq0,genesr0,geness0, genest0),
#                  cbind(genesq1,genesr1,geness1, genest1))  
#   
#   colnames(genes) <- paste0("Gene", 1:p)
#   rownames(genes) <- paste0("Subject",1:n)
#   
#   # not used for now
#   # The coefficients are constructed to have alternating signs and to be exponentially
#   # decreasing. (Friedman et al 2010, Journal of Statistical Software)
#   # beta <- (-1)^{1:p} * exp(-(2*1:p-1)/100)
#   # beta <- (-1)^{1:p} * exp(-(2*1:p - 1 )/600)*sin(1:p)
#   
#   y.star <- genes %*% beta.genes
#   #as.matrix(X) %>% dim
#   #dim(genes)
#   
#   error <- rnorm(n)
#   k <- sqrt(var(y.star)/(signal.to.noise.ratio*var(error)))
#   
#   y <- y.star + k*error 
#   
#   result <- as.matrix(cbind(y,E,genes))
#   class(result) <- append(class(result), "eset")
#   
#   return(result)
# }



## ---- heatmap-expression ----
plot.eset <- function(x, color.palette = colorRampPalette(rev(brewer.pal(n = 7, name = "Reds")))(100),
                      n = 400, n0 = 200, xlab = "Subjects", ylab = "Gene Expression", ...){
  
  
  # function to generate heatmap with x and y labels
  # note that you need to use this hack because pheatmap doesnt have 
  # a default for X and Y labels
  # data: matrix of gene expression data (rows are subjects, columns are genes)
  # takes object of class expression
  # col.pal: RColorBrewer color palette  
  genes <- grep("Gene(\\d+)$",x  %>% colnames(), perl = T, value = T)
  annotation_col = data.frame(Exposure = factor(c(rep("E=0",n0),rep("E=1", n-n0))))
  rownames(annotation_col) = paste0("Subject", 1:n)
    
  setHook("grid.newpage", 
          function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), 
          action="prepend")
  pheatmap::pheatmap(t(x[,genes]) %>% magrittr::set_colnames(paste0("Subject", 1:n)), 
                     annotation_col = annotation_col,
                     color = color.palette,
                     show_colnames = F,
                     show_rownames = F,
                     cluster_cols = F, 
                     cluster_rows = F, 
                     gaps_col = n0, ...)
  setHook("grid.newpage", NULL, "replace")
  grid.text(xlab, y=-0.07, gp=gpar(fontsize=16))
  grid.text(ylab, x=-0.07, rot=90, gp=gpar(fontsize=16))
}


## ---- heatmap-correlation ----

plot.similarity <- function(x, 
                          color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(100), ...){
  
  # function to generate heatmap
  # data: matrix of true correlation (P x P matrix where P is the number of genes)
  # takes object of class truecorr
  # col.pal: RColorBrewer color palette  

  pheatmap::pheatmap(x, 
                     color = color,
                     show_colnames = F,
                     show_rownames = F,
                     cluster_cols = F, 
                     cluster_rows = F,
#                      breaks = seq(min(min_max_heat$V2), max(min_max_heat$V1), length.out = 101) ,
#                      legend_breaks = round(seq(min(min_max_heat$V2), max(min_max_heat$V1), length.out = 12),1),
#                      legend_labels = round(seq(min(min_max_heat$V2), max(min_max_heat$V1), length.out = 12),1), 
                     drop_levels = FALSE, ...)
}


## ---- new-pc ----

new.pc <- function(data,loadings) {
  as.matrix(data) %*% as.matrix(loadings)
}


## ---- summary-se ----

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     lower = quantile(xx[[col]], 0.025,na.rm=na.rm),
                     upper = quantile(xx[[col]],0.975,na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  #datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}


## ---- multiplot ----

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



## ---- stability ----
my_cl_valid <- function(rowDel, data = X.train, genes.x0 = genes.x0, genes.x1 = genes.x1,
                        clust.dist.type = args[1], cluster = clusters$cluster,
                        distance = distance){
  
  #   data = X.train; genes.x0 = genes.x0; genes.x1 = genes.x1;
  #   clust.dist.type = args[1]; rowDel=1 ; cluster = clusters$cluster
  # the clValid::stability function takes rows as the items to be clustered
  # therefore rows are genes and columns are subjects
  
  
  genes.x0.del <- genes.x0[-rowDel,]
  genes.x1.del <- genes.x1[-rowDel,]
  
  distanceDel <- switch(clust.dist.type,
                        corr = {WGCNA::cor(data[-rowDel,])},
                        corr0 = {WGCNA::cor(genes.x0.del)},
                        corr1 = {WGCNA::cor(genes.x1.del)},
                        diffcorr = {abs(WGCNA::cor(genes.x1.del)-WGCNA::cor(genes.x0.del))},
                        difftom = {abs(WGCNA::TOMsimilarityFromExpr(genes.x1.del, corType = "pearson", power = 6)-
                                         WGCNA::TOMsimilarityFromExpr(genes.x0.del, corType = "pearson", power = 6))},
                        tom0 = {WGCNA::TOMsimilarityFromExpr(genes.x0.del, corType = "pearson", power = 6)},
                        tom1 = {WGCNA::TOMsimilarityFromExpr(genes.x1.del, corType = "pearson", power = 6)},
                        tom = {WGCNA::TOMsimilarityFromExpr(data[-rowDel,],  corType = "pearson", power = 6)})
  
  clDel <- hclust(as.dist(if (clust.dist.type %in% c("diffcorr","difftom")) distanceDel else 1-distanceDel), method = "average")
  cuttreeDel <- dynamicTreeCut::cutreeDynamic(clDel, 
                                              distM = if (clust.dist.type %in% c("diffcorr","difftom")) distanceDel else 1-distanceDel, 
                                              deepSplit = 1)
  
  # check if all cluster groups are 0 which means no cluster assignment and everyone is in their own group
  clustersDel <- data.table(gene = paste0("Gene",1:p),cluster = if (all(cuttreeDel==0)) 1:p else cuttreeDel)
  setkey(clustersDel,"cluster")
  
  stab <- clValid::stability(mat = t(X.train),
                             Dist = as.dist(if (args[1] %in% c("diffcorr","difftom")) distance else 1-distance),
                             del = rowDel,
                             cluster = cluster, 
                             clusterDel = clustersDel$cluster)
  
  return(stab)
}


## ---- Analysis-Functions ----

# filtering on univariate p-value, taking the lowet XX percent and then running 
# multiple linear regression on it
uni_fun <- function(train, 
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
    
    # False Positive Rate
    uni.FPR <- sum(c(gene.names,uni.S.hat) %ni% s0)/(p - length(s0))
    
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
                      k_means = NULL){
  
  #   result %>% names
#               stability = F; gene_groups = result[["clusters"]]; 
#               x_train = result[["X_train"]] ; x_test = result[["X_test"]] ; y_train = result[["Y_train"]] ; y_test = result[["Y_test"]]; 
#               filter = F; filter_var = F; include_E = F; include_interaction = F; 
#               s0 = result[["S0"]]; p = 1000 ;
#               model = "lasso"; summary = "avg"; topgenes = NULL;
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
  
  print(paste(summary,model,"filter = ", filter, "filter_var = ",filter_var, "include_E = ", include_E, "include_interaction = ", include_interaction, sep = ","))
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
  
  gene.names <- colnames(x_train_mod)[which(colnames(x_train_mod) %ni% "E")]
  
  # the number of clusters in the modified set
  # remove genes which have a cluster assignment of 0, meaning they dont cluster well
  n.clusters <-   gene_groups[cluster != 0][gene %in% gene.names]$cluster %>% unique %>% length
  cluster_names <- gene_groups[cluster != 0][gene %in% gene.names]$cluster %>% unique
  
  
  clustered.genes <- lapply(cluster_names, function(i) {
    x_train_mod[,intersect(gene_groups[cluster == i]$gene, gene.names), drop = FALSE] %>% scale(center = T, scale = F) %>% as.matrix
  }
  )
  
  # remove clusters that have no data in them due to filtering
  clustered.genes <- Filter(function(i) ncol(i) > 0, clustered.genes)
  print(sapply(clustered.genes, dim))
  
  clust_data <- switch(summary, 
                       avg = plyr::ldply(clustered.genes, rowMeans) %>% t, 
                       pc = {
                         # output from prcomp
                         pc.train <- lapply(clustered.genes, function(i) prcomp(i, center = F))
                         # extract 1st PC for each cluster
                         plyr::ldply(pc.train, function(i) i$x[,1]) %>% t
                       },
                       spc = {
                         out <- lapply(clustered.genes, function(i) PMA::SPC.cv(i, 
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
  
  if (summary == "spc") {
    out_best <- clust_data[[2]] 
    out_best_v <- clust_data[[3]] # used to calculate S.hat (number of non zero components) 
    
    # these are the gene names corresponding to the non zero loadings from sparsePCA
    non_zero_pc_components <- lapply(out_best_v, function(i) 
      i[which(i != 0), ,drop = F] %>% rownames) %>% unlist
  }
  
  clust_data <- if (summary == "spc") clust_data[[1]] else clust_data
  colnames(clust_data) <- paste0(summary,cluster_names)
  
  #head(clust_data)  
  
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
  
  X.model.formula <- model.matrix(model.formula, data = if (include_E) cbind(clust_data,x_train_mod[,"E", drop = F]) else clust_data %>% as.data.frame )
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
                  }
  )
  coefs
  if (stability) {
    # remove intercept for stability measures
    return(coefs %>% magrittr::extract(-1, , drop = F))
  } else {
    
    non_zero_clusters <- coefs[-1, , ][coef.est != 0] %>% 
      magrittr::use_series("Gene")
    
    n.non_zero_clusters <- coefs[-1, , ][coef.est != 0] %>% 
      magrittr::use_series("Gene") %>% 
      length
    
    # test data
    x_test_mod = if (filter & !include_E) {x_test[, topgenes] %>% as.data.frame} else if (!filter & include_E) {
      x_test %>% as.data.frame } else if (!filter & !include_E) {
        x_test[,which(colnames(x_test) %ni% "E")] %>% as.data.frame} else if (filter & include_E) {
          x_test[, c(topgenes,"E")] %>% as.data.frame
        }
    
    clustered.genes_test <- lapply(cluster_names, function(i) {
      x_test_mod[,intersect(gene_groups[cluster == i]$gene, gene.names), drop = FALSE] %>% scale(center = T, scale = F) %>% as.matrix
    }
    )
    
    # remove clusters that have no data in them due to filtering
    clustered.genes_test <- Filter(function(i) ncol(i) > 0, clustered.genes_test)
    
    clust_data_test <- switch(summary, 
                              avg = plyr::ldply(clustered.genes_test, rowMeans) %>% t, 
                              pc = {
                                # output from prcomp
                                pc.train <- lapply(clustered.genes, function(i) prcomp(i, center = F))
                                # extract loadings for 1st PC i.e. 1st eigenvector of X, for each cluster
                                V <- lapply(pc.train, function(i) i$rotation[,1])
                                mapply(new.pc, data = clustered.genes_test, 
                                       loadings = V)
                              },
                              spc = {
                                V <- lapply(out_best, function(i) i$v)
                                mapply(new.pc, data = clustered.genes_test,
                                       loadings = V)
                              }
    )
    
    colnames(clust_data_test) <- paste0(summary,cluster_names)
    
    # need intercept for prediction
    model.formula_test <- if (include_interaction & include_E) {
      paste0("~ 1+(",paste0(colnames(clust_data_test), collapse = "+"),")*E") %>% as.formula} else if (!include_interaction & include_E) {
        paste0("~1+",paste0(colnames(clust_data_test), collapse = "+"),"+E") %>% as.formula} else if (!include_interaction & !include_E) {
          paste0("~1+",paste0(colnames(clust_data_test), collapse = "+")) %>% as.formula
        }
    
    X.model.formula_test <- model.matrix(model.formula_test, 
                                         data = if (include_E) cbind(clust_data_test,x_test_mod[,"E", drop = F]) else clust_data_test %>% as.data.frame)
    
    # need to get the genes corresponding to the non-zero clusters
    # NOTE: this also includes non-zero cluster:Environment interactions
    
    # genes corresponding to non-zero clusters
    # return non zero sparse loadings if spc combined with lm, else
    # return non-zero clusters for penalization methods, even if your using spc
    clust.S.hat <- if (summary == "spc" & model == "lm") non_zero_pc_components else {gene_groups %>%  
                                                                                        magrittr::extract(cluster %in% (grep("[0-9]",non_zero_clusters, value = T) %>% 
                                                                                                                          gsub(summary,"", .) %>% 
                                                                                                                          gsub(":E","",.) %>%
                                                                                                                          unique %>% 
                                                                                                                          as.numeric)) %>% magrittr::use_series("gene") }
    
    
    
    # True Positive Rate
    clust.TPR <- length(intersect(clust.S.hat, s0))/length(s0)
    
    # False Positive Rate
    clust.FPR <- sum(clust.S.hat %ni% s0)/(p - length(s0))
    
    # Mean Squared Error
    clust.mse <- crossprod(X.model.formula_test %*% coefs$coef.est - y_test)/length(y_test)  
    
    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)
    
    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    clust.r2 <- (mse_null - clust.mse)/mse_null
    
    clust.adj.r2 <- 1 - (1 - clust.r2)*(nrow(x_test) - 1)/(nrow(x_test) - n.non_zero_clusters - 1)
    
    
    ls <- list(clust.mse = clust.mse, clust.r2 = clust.r2, clust.adj.r2 = clust.adj.r2, clust.S.hat = length(clust.S.hat), 
               clust.TPR = clust.TPR, clust.FPR = clust.FPR)
    
    names(ls) <- c(paste0(if (is.null(k_means)) "clust_" else "Eclust_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0(if (is.null(k_means)) "clust_" else "Eclust_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_r2"),
                   paste0(if (is.null(k_means)) "clust_" else "Eclust_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_adjr2"),
                   paste0(if (is.null(k_means)) "clust_" else "Eclust_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_Shat"),
                   paste0(if (is.null(k_means)) "clust_" else "Eclust_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_TPR"),
                   paste0(if (is.null(k_means)) "clust_" else "Eclust_",summary,"_",model,ifelse(include_interaction,"_yes","_no"),"_FPR"))
    return(ls)
    
  }
}


pen_fun <- function(x_train, 
                    x_test, 
                    y_train, 
                    y_test, 
                    s0,
                    model,
                    true_beta,
                    topgenes = NULL,
                    stability = F,
                    filter = F, 
                    include_E = F, 
                    include_interaction = F,
                    p = 1000,
                    filter_var = F){
  
  #   stability = F; x_train = result[["X_train"]] ; x_test = result[["X_test"]] ; 
  #   y_train = result[["Y_train"]] ; y_test = result[["Y_test"]]; 
  #   filter = F; filter_var = F; include_E = F; include_interaction = F; 
  #   s0 = result[["S0"]]; p = 1000 ;
  #   model = "scad"; topgenes = NULL; true_beta = result[["beta_truth"]]
  
  #     stability = F; x_train = result_interaction[["X_train"]] ; x_test = result_interaction[["X_test"]] ; 
  #     y_train = result_interaction[["Y_train"]] ; y_test = result_interaction[["Y_test"]]; 
  #     filter = F; filter_var = F; include_E = T; include_interaction = T; 
  #     s0 = result_interaction[["S0"]]; p = 1000 ;
  #     model = "scad"; topgenes = NULL; true_beta = result_interaction[["beta_truth"]]
  
  # model: "scad", "mcp", "lasso", "elasticnet", "ridge" 
  # filter: T or F based on univariate filter
  
  print(paste(model,"filter = ", filter, "filter_var = ",filter_var, "include_E = ", include_E, "include_interaction = ", include_interaction, sep = ","))
  
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
  
  #gene.names <- colnames(x_train)[which(colnames(x_train) %ni% "E")]
  
  # penalization model
  pen_model <- switch(model, 
                      lasso = glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 1),
                      elasticnet = glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 0.5),
                      ridge = glmnet::cv.glmnet(x = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train, alpha = 0),
                      scad = ncvreg::cv.ncvreg(X = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train,
                        family = "gaussian", penalty = "SCAD"),
                      mcp = ncvreg::cv.ncvreg(X = if (!include_E) as.matrix(x_train[,-grep("E", colnames(x_train))]) else
                        as.matrix(x_train), y = y_train,
                        family = "gaussian", penalty = "MCP")
  )
  
  # oracle penalization model
  pen_model_oracle <- switch(model, 
                             lasso = glmnet::cv.glmnet(x = as.matrix(x_train[,s0]), y = y_train, alpha = 1),
                             elasticnet = glmnet::cv.glmnet(x = as.matrix(x_train[,s0]), y = y_train, alpha = 0.5),
                             ridge = glmnet::cv.glmnet(x = as.matrix(x_train[,s0]), y = y_train, alpha = 0),
                             scad = ncvreg::cv.ncvreg(X = x_train[,s0], y = y_train,
                                                      family = "gaussian", penalty = "SCAD"),
                             mcp = ncvreg::cv.ncvreg(X = x_train[,s0], y = y_train,
                                                     family = "gaussian", penalty = "MCP")
  )
  
  # here we give the coefficient stability on the individual genes
  coefs <- coef(pen_model, s = "lambda.min") %>% 
    as.matrix %>% 
    as.data.table(keep.rownames = TRUE) %>% 
    magrittr::set_colnames(c("Gene","coef.est")) %>% 
    magrittr::extract(-1,)
  
  
  if (stability) {
    # remove intercept for stability measures
    return(coefs)
  } else {
    
    pen.S.hat <- coefs[coef.est != 0] %>% magrittr::use_series("Gene") 
    
    pen.pred <- if (model %in% c("lasso","elasticnet","ridge")) {
      predict(pen_model, newx =  if (!include_E) as.matrix(x_test[,-grep("E", colnames(x_test))]) else
        as.matrix(x_test), s = "lambda.min") } else if (model %in% c("scad","mcp")) { 
          predict(pen_model, X =  if (!include_E) as.matrix(x_test[,-grep("E", colnames(x_test))]) else
            as.matrix(x_test), 
            lambda = pen_model$lambda.min)
        }
    
    pen.pred.oracle <- if (model %in% c("lasso","elasticnet","ridge")) {
      predict(pen_model_oracle, newx = x_test[,s0], s = "lambda.min") } else if (model %in% c("scad","mcp")) { 
        predict(pen_model_oracle, X = x_test[,s0], 
                lambda = pen_model_oracle$lambda.min)
      }
    
    # Mean Squared Error
    pen.mse <- crossprod(pen.pred - y_test)/length(y_test)
    
    # Mean Squared Error Oracle
    pen.mse.oracle <- crossprod(pen.pred.oracle - y_test)/length(y_test)
    
    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)
    
    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    pen.r2 <- (mse_null - pen.mse)/mse_null
    
    pen.adj.r2 <- 1 - (1 - pen.r2)*(length(y_test) - 1)/(length(y_test) - length(pen.S.hat) - 1)
    
    # True Positive Rate
    pen.TPR <- length(intersect(pen.S.hat, s0))/length(s0)
    
    # False Positive Rate
    pen.FPR <- sum(pen.S.hat %ni% s0)/(p - length(s0))
    
    # model error
    identical(true_beta %>% rownames(),coefs[["Gene"]])
    pen.model.error <- {(true_beta - coefs[["coef.est"]]) %>% t} %*% WGCNA::cor(x_test[,coefs[["Gene"]]]) %*% (true_beta - coefs[["coef.est"]])
    
    ls <- list(pen.mse = as.numeric(pen.mse), pen.r2 = as.numeric(pen.r2), 
               pen.adj.r2 = as.numeric(pen.adj.r2), pen.S.hat = length(pen.S.hat), 
               pen.TPR = pen.TPR, pen.FPR = pen.FPR, pen.relative.mse = pen.mse/pen.mse.oracle, pen.model.error = pen.model.error)
    names(ls) <- c(paste0("pen","_",model,ifelse(filter_var,"_na","_na"),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0("pen","_",model,ifelse(filter_var,"_na","_na"),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_r2"),
                   paste0("pen","_",model,ifelse(filter_var,"_na","_na"),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_adjr2"),
                   paste0("pen","_",model,ifelse(filter_var,"_na","_na"),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_Shat"),
                   paste0("pen","_",model,ifelse(filter_var,"_na","_na"),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_TPR"),
                   paste0("pen","_",model,ifelse(filter_var,"_na","_na"),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_FPR"),
                   paste0("pen","_",model,ifelse(filter_var,"_na","_na"),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_relmse"),
                   paste0("pen","_",model,ifelse(filter_var,"_na","_na"),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_modelerror"))
    return(ls)
  }
  
}

group_pen_fun <- function(x_train, 
                          x_test, 
                          y_train, 
                          y_test, 
                          s0,
                          model,
                          true_beta,
                          gene_groups,
                          topgenes = NULL,
                          stability = F,
                          filter = F, 
                          include_E = F, 
                          include_interaction = F,
                          p = 1000,
                          filter_var = F){
  
  #   stability = F; x_train = result[["X_train"]] ; x_test = result[["X_test"]] ; 
  #   y_train = result[["Y_train"]] ; y_test = result[["Y_test"]]; 
  #   filter = F; filter_var = F; include_E = F; include_interaction = F; 
  #   s0 = result[["S0"]]; p = 1000 ;
  #   model = "gglasso"; topgenes = NULL; true_beta = result[["beta_truth"]]
  #   gene_groups = result[["clusters"]]
  #   
  #   # interaction
  #   stability = F; x_train = result_interaction[["X_train"]] ; x_test = result_interaction[["X_test"]] ; 
  #   y_train = result_interaction[["Y_train"]] ; y_test = result_interaction[["Y_test"]]; 
  #   filter = F; filter_var = F; include_E = T; include_interaction = T; 
  #   s0 = result_interaction[["S0"]]; p = 1000 ;
  #   model = "gglasso"; topgenes = NULL; true_beta = result_interaction[["beta_truth"]]
  #   gene_groups = result_interaction[["gene_groups_inter"]]
  
  # model: "gglasso"
  # filter: T or F based on univariate filter
  
  print(paste(model,"filter = ", filter, "filter_var = ",filter_var, "include_E = ", include_E, "include_interaction = ", include_interaction, sep = ","))
  
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
  
  gene.names <- colnames(x_train)[which(colnames(x_train) %ni% "E")]
  
  #   cv.glasso <- gglasso::cv.gglasso(as.matrix(x_train[,gene_groups[cluster != 0]$gene]), y_train, group = gene_groups[cluster != 0]$cluster,
  #                           loss = "ls")
  #   
  #   coef.glasso <- coef(cv.glasso, s = "lambda.min") %>% as.data.table(keep.rownames = T)
  #   setnames(coef.glasso, c("gene","coef"))
  #   fit.glasso <- predict(cv.glasso, newx = as.matrix(DT.test[,clusters[cluster!=0]$gene]), s = "lambda.min")
  #   
  #   # Mean Squared Error
  #   gglasso.mse <- crossprod(fit.glasso - DT.test[,"V1"])/nrow(DT.test)
  #   
  #   # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
  #   gglasso.r2 <- (mse.null-gglasso.mse)/mse.null
  
  
  
  # penalization model
  grp_pen_model <- switch(model, 
                          gglasso =   if (include_interaction) {
                            gglasso::cv.gglasso(as.matrix(x_train[,gene_groups[cluster != 0]$gene]), 
                                                y_train,
                                                pf = gene_groups %>% distinct(cluster) %>% magrittr::use_series("pf"),
                                                group = gene_groups[cluster != 0]$cluster,
                                                loss = "ls")} else {
                                                  gglasso::cv.gglasso(as.matrix(x_train[,gene_groups[cluster != 0]$gene]), 
                                                                      y_train,
                                                                      group = gene_groups[cluster != 0]$cluster,
                                                                      loss = "ls")
                                                }
  )
  
  # oracle penalization model
  grp_pen_model_oracle <- switch(model, 
                                 gglasso =  gglasso::cv.gglasso(as.matrix(x_train[,gene_groups[cluster != 0][gene %in% s0]$gene]), 
                                                                y_train, 
                                                                group = factor(gene_groups[cluster != 0][gene %in% s0]$cluster) %>% as.numeric,
                                                                loss = "ls"))
  
  # here we give the coefficient stability on the individual genes
  coefs <- coef(grp_pen_model, s = "lambda.min") %>% 
    as.matrix %>% 
    as.data.table(keep.rownames = TRUE) %>% 
    magrittr::set_colnames(c("Gene","coef.est")) %>% 
    magrittr::extract(-1,)
  
  if (stability) {
    return(coefs)
  } else {
    
    grp.pen.S.hat <- coefs[Gene != "E"][coef.est != 0] %>% magrittr::use_series("Gene") %>% 
      gsub(":E","",.) %>% unique
    
    grp.pen.pred <- predict(grp_pen_model, newx = as.matrix(x_test[,gene_groups[cluster !=0 ]$gene]), s = "lambda.min")
    grp.pen.pred.oracle <- predict(grp_pen_model_oracle, newx = as.matrix(x_test[,s0]), s = "lambda.min")
    # this is to make sure that the order of coefficients corresponds to that of the test data:
    # identical(coef(grp_pen_model)[-1,,drop=F] %>% rownames(),gene_groups[cluster !=0 ]$gene )
    
    # Mean Squared Error
    grp.pen.mse <- crossprod(grp.pen.pred - y_test)/length(y_test)
    
    # Mean Squared Error Oracle
    grp.pen.mse.oracle <- crossprod(grp.pen.pred.oracle - y_test)/length(y_test)
    
    # mse.null
    mse_null <- crossprod(mean(y_test) - y_test)/length(y_test)
    
    # the proportional decrease in model error or R^2 for each scenario (pg. 346 ESLv10)
    grp.pen.r2 <- (mse_null - grp.pen.mse)/mse_null
    
    grp.pen.adj.r2 <- 1 - (1 - grp.pen.r2)*(length(y_test) - 1)/(length(y_test) - length(grp.pen.S.hat) - 1)
    
    # True Positive Rate
    grp.pen.TPR <- length(intersect(grp.pen.S.hat, s0))/length(s0)
    
    # False Positive Rate
    grp.pen.FPR <- sum(grp.pen.S.hat %ni% s0)/(p - length(s0))
    
    # model error
    identical(true_beta %>% rownames(),coefs[["Gene"]])
    
    # this is required to be able to take proper differences between vectors, so we join them first
    truth <- true_beta %>% as.data.table(keep.rownames = T) %>% magrittr::set_colnames(c("Gene","true_beta")) %>% setkey(Gene)
    setkey(coefs,Gene)
    
    tmp <- truth[coefs]
    tmp[,diff:=true_beta-coef.est]
    
    grp.pen.model.error <- {(tmp$diff) %>% t} %*% WGCNA::cor(x_test[,tmp[["Gene"]]]) %*% (tmp$diff)
    
    ls <- list(grp.pen.mse = as.numeric(grp.pen.mse), grp.pen.r2 = as.numeric(grp.pen.r2), 
               grp.pen.adj.r2 = as.numeric(grp.pen.adj.r2), grp.pen.S.hat = length(grp.pen.S.hat), 
               grp.pen.TPR = grp.pen.TPR,
               grp.pen.FPR = grp.pen.FPR,
               grp.pen.relative.mse = grp.pen.mse/grp.pen.mse.oracle, 
               #grp.pen.relative.mse = NA,
               grp.pen.model.error = grp.pen.model.error)
    names(ls) <- c(paste0("group_na","_",model,ifelse(filter_var,"",""),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_mse"),
                   paste0("group_na","_",model,ifelse(filter_var,"",""),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_r2"),
                   paste0("group_na","_",model,ifelse(filter_var,"",""),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_adjr2"),
                   paste0("group_na","_",model,ifelse(filter_var,"",""),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_Shat"),
                   paste0("group_na","_",model,ifelse(filter_var,"",""),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_TPR"),
                   paste0("group_na","_",model,ifelse(filter_var,"",""),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_FPR"),
                   paste0("group_na","_",model,ifelse(filter_var,"",""),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_relmse"),
                   paste0("group_na","_",model,ifelse(filter_var,"",""),ifelse(include_E,"",""),ifelse(include_interaction,"_yes","_no"),"_modelerror"))
    return(ls)
  }
  
}

## ---- test-pos-def ----

posdef <- function(x){
if (all(eigen(x)$values>0)) paste("Matrix is positive definite") else paste("Matrix is not pos-def")
}


"%ni%" <- Negate("%in%")


## ---- proftable ----

#' An alternative to \code{summaryRprof()}
#' 
#' \code{proftools} parses a profiling file and prints an easy-to-understand
#' table showing the most time-intensive function calls. 
#' 
#' Line numbers are included if \code{Rprof()} was run with 
#' \code{line.numbering=TRUE}. If it was run with \code{memory.profiling=TRUE},
#' this function will probably break.
#' 
#' Below the table are printed any files identified if line numbering is true,
#' the total time recorded by \code{Rprof()}, and the "parent call".  The
#' parent call consists of the parent call stack of all the call stacks in the\
#' table. Note that this is the parent call stack of only the printed lines,
#' not of all stacks recorded by \code{Rprof()}. This makes the table easier to read and fit into the console. 
#' 
#' @export
#' @param file A profiling file generated by \code{Rprof()}
#' @param lines The number of lines (call stacks) you want returned. Lines are
#' printed from most time-intensive to least.
proftable <- function(file, lines = 10) {
  profdata <- readLines(file)
  interval <- as.numeric(strsplit(profdata[1L], "=")[[1L]][2L]) / 1e+06
  filelines <- grep("#File", profdata)
  files <- profdata[filelines]
  profdata <- profdata[-c(1, filelines)]
  total.time <- interval * length(profdata)
  ncalls <- length(profdata)
  profdata <- gsub("\\\"| $", "", profdata)
  calls <- lapply(profdata, function(x) rev(unlist(strsplit(x, " "))))
  stacktable <- as.data.frame(table(sapply(calls, function(x) paste(x, collapse = " > "))) / ncalls * 100, stringsAsFactors = FALSE)
  stacktable <- stacktable[order(stacktable$Freq[], decreasing = TRUE), 2:1]
  colnames(stacktable) <- c("PctTime", "Call")
  stacktable <- head(stacktable, lines)
  shortcalls = strsplit(stacktable$Call, " > ")
  shortcalls.len <- range(sapply(shortcalls, length))
  parent.call <- unlist(lapply(seq(shortcalls.len[1]), function(i) Reduce(intersect, lapply(shortcalls,"[[", i))))
  shortcalls <- lapply(shortcalls, function(x) setdiff(x, parent.call))
  stacktable$Call = sapply(shortcalls, function(x) paste(x, collapse = " > "))
  if (length(parent.call) > 0) {
    parent.call <- paste(paste(parent.call, collapse = " > "), "> ...")
  } else {
    parent.call <- "None"
  }
  frac <- sum(stacktable$PctTime)
  attr(stacktable, "total.time") <- total.time
  attr(stacktable, "parent.call") <- parent.call
  attr(stacktable, "files") <- files
  attr(stacktable, "total.pct.time") <- frac
  print(stacktable, row.names=FALSE, right=FALSE, digits=3)
  if(length(files) > 0) {
    cat("\n")
    cat(paste(files, collapse="\n"))
    cat("\n")
  }
  cat(paste("\nParent Call:", parent.call))
  cat(paste("\n\nTotal Time:", total.time, "seconds\n"))
  cat(paste0("Percent of run time represented: ", format(frac, digits=3)), "%")
  
  invisible(stacktable)
}
