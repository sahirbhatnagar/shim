##################################
# R source code file for conducting simulations
# on HYDRA cluster
# Git: this is on the eclust repo, simulation branch
# Created by Sahir,  April 2, 2016
# Updated:
# Notes:
# This is a modified simulation and different from simulation1
# Its based on code from Network analysis book by Horvath
##################################

rm(list=ls())
source("packages.R")
source("functions.R")

parametersDf <- expand.grid(rho = c(0.1,0.35,0.75,0.95),
                            p = c(500, 1000, 3000),
                            n = 500, n0 = 250,
                            cluster_distance = c("corr"),
                            rhoOther = 0.6,
                            betaMean = 4,
                            betaE = 5,
                            alphaMean = 2,
                            nActive = 50,
                            includeInteraction = TRUE,
                            includeStability = TRUE,
                            distanceMethod = "euclidean",
                            clustMethod = "hclust",
                            cutMethod = "dynamic",
                            method = "complete",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

#str(parametersDf)
#parametersDf <- transform(parametersDf, nBlocks = p/blocksize)
parameterIndex <- commandArgs(trailingOnly = T)

parameterIndex = 4
simulationParameters <- parametersDf[parameterIndex,, drop = F]

## ---- generate-data ----
message("generating data")

p <- simulationParameters[,"p"];
n <- simulationParameters[,"n"];
n0 <- simulationParameters[,"n0"];
n1 <- n - n0
cluster_distance <- simulationParameters[,"cluster_distance"]
rhoOther <- simulationParameters[,"rhoOther"];
betaMean <- simulationParameters[,"betaMean"];
betaE <- simulationParameters[,"betaE"]
alphaMean <- simulationParameters[,"alphaMean"];
nBlocks <- simulationParameters[,"nBlocks"];
rho <- simulationParameters[,"rho"];
nActive <- simulationParameters[,"nActive"];
includeInteraction <- simulationParameters[,"includeInteraction"]
includeStability <- simulationParameters[,"includeStability"]
distanceMethod <- simulationParameters[,"distanceMethod"]
clustMethod <- simulationParameters[,"clustMethod"]
cutMethod <- simulationParameters[,"cutMethod"]
method <- simulationParameters[,"method"]
K.max <- simulationParameters[,"K.max"]
B <- simulationParameters[,"B"]

# in this simulation its blocks 3 and 4 that are important
d0 <- simModule(n = n0, p = p, rho = c(rho,-rho), exposed = FALSE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.3,
                maxCor = 1,
                corPower = 0.3,
                propNegativeCor = 0.1,
                backgroundNoise = 0.2,
                signed = TRUE,
                leaveOut = 3:4)

d1 <- simModule(n = n1, p = p, rho = c(rho, rho), exposed = TRUE,
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
t1 <- table(truemodule1)
table(truemodule0,truemodule1)

# Convert labels to colors for plotting
moduleColors <- labels2colors(truemodule1)

X <- rbind(d0$datExpr, d1$datExpr) %>%
  magrittr::set_colnames(paste0("Gene", 1:p)) %>%
  magrittr::set_rownames(paste0("Subject",1:n))

dim(X)


betaMainEffect <- vector("double", length = p)
betaMainInteractions <- vector("double", length = p)
# first assign random uniform to every gene in cluster 3 and 4,
# then randomly remove so that thers only nActive left
betaMainEffect[which(truemodule1 %in% 3:4)] <- runif(sum(truemodule1 %in% 3:4),
                                                     betaMean - 0.1, betaMean + 0.1)

# randomly set some coefficients to 0 so that there are only nActive non zero
betaMainEffect <- replace(betaMainEffect,
                          sample(which(truemodule1 %in% 3:4), sum(truemodule1 %in% 3:4) - nActive,
                                 replace = FALSE), 0)

betaMainInteractions[which(betaMainEffect!=0)] <- runif(nActive, alphaMean - 0.1, alphaMean + 0.1)

beta <- c(betaMainEffect,
          betaE,
          betaMainInteractions)
plot(beta)

result <- generate_data(p = p, n = n, n0 = n0, X = X,
                        beta = beta, include_interaction = includeInteraction,
                        cluster_distance = cluster_distance,
                        EclustAddDistance = "fisherScore",
                        signal_to_noise_ratio = 1/2,
                        distanceMethod = distanceMethod,
                        clustMethod = clustMethod,
                        cutMethod = cutMethod,
                        method = method,
                        K.max = K.max, B = B)


## ---- univariate-pvalue -----

message("Starting univariate p-value with interaction")

# filtering on univariate p-value, taking the lowet 5 percent and then running
# multiple linear regression on it
# the correlations are really messing up the oracle model which is why when
# you filter you get better fit
# stability = TRUE makes the function return coefficients only. use only for
# calculating stability measures
# output is in the following format: eg. uni_na_lm_yes_mse
# which represents method_clustermeasure_model_includeE_measure
uni_res <- uniFit(train = result[["DT_train"]], test = result[["DT_test"]],
                  percent = 0.05, stability = F,
                  include_E = T,
                  include_interaction = includeInteraction,
                  filter_var = F, p = p,
                  s0 = result[["S0"]], true_beta = result[["beta_truth"]])

if (includeStability) {
  uni_stab <- mapply(uniFit,
                     train = result[["DT_train_folds"]],
                     MoreArgs = list(test = result[["DT_test"]],
                                     s0 = result[["S0"]],
                                     true_beta = result[["beta_truth"]],
                                     stability = T,
                                     include_E = T,
                                     include_interaction = includeInteraction,
                                     percent = 0.05,
                                     filter_var = F,
                                     p = p),
                     SIMPLIFY = F)

  # Make the combinations of list elements
  ll <- combn(uni_stab, 2, simplify = F)

  # Pairwise correlations of the model coefficients for each of the 10 CV folds
  uni_mean_stab <- lapply(c("pearson","spearman"), function(i) {
    res <- sapply(ll,
                  function(x) WGCNA::cor(x[[1]]$coef.est,
                                         x[[2]]$coef.est,
                                         method = i,use = 'pairwise.complete.obs')
    ) %>% mean
    names(res) <- paste0("uni_na_lm_yes_",i)
    return(res)
  }
  )

  # Jaccard index
  uni_jacc <- sapply(ll, function(x) {
    A = x[[1]][coef.est != 0]$Gene
    B = x[[2]][coef.est != 0]$Gene
    length(intersect(A,B))/length(union(A,B))
  }
  ) %>% mean

}

message("done univariate p-value with interaction")








## ---- cluster-and-regress----

# we will treat the clusters as fixed i.e., even if we filter, or
# do cross validation, the group labels are predetermined by the
# above clustering procedure

print("starting cluster and regress with interaction")

clust_res <- mapply(clust_fun,
                    #summary = rep(c("pc","spc","avg"), each = 3),
                    #model = rep(c("lm", "lasso","elasticnet"), 3),
                    summary = c("avg","avg"),
                    model = c("lasso","shim"),
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
                                    gene_groups = result[["clustersAll"]],
                                    clust_type = "clust"),
                    SIMPLIFY = F,
                    USE.NAMES = F)

result %>% names
clust_res %>% unlist

if (includeStability) {
  clust_stab <- mapply(function(summary,
                                model) mapply(clust_fun,
                                              x_train = result[["X_train_folds"]],
                                              y_train = result[["Y_train_folds"]],
                                              MoreArgs = list(stability = T,
                                                              summary = summary,
                                                              model = model,
                                                              filter = F,
                                                              filter_var = F,
                                                              include_E = T,
                                                              include_interaction = includeInteraction,
                                                              gene_groups = result[["clusters"]],
                                                              p = p),
                                              SIMPLIFY = F),
                       #summary = rep(c("pc","spc","avg"), each = 3),
                       #model = rep(c("lm", "lasso","elasticnet"), 3),
                       summary =  c("avg","avg"),
                       model = c("lasso","shim"),
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
                       summary = c("avg","avg"),
                       model = c("lasso","shim"),
                       USE.NAMES = F)


  # Pairwise correlations of the model coefficients for each of the 10 CV folds
  clust_mean_stab <- lapply(seq_along(ll), function(j) {
    lapply(c("pearson","spearman"), function(i) {
      res <- sapply(ll[[j]] , function(x) WGCNA::cor(x[[1]]$coef.est, x[[2]]$coef.est, method = i)) %>% mean
      names(res) <- paste0(clust_labs[[j]], i)
      return(res)
    }
    )
  }
  )

  clust_mean_stab %>% unlist

  # Jaccard index
  clust_jacc <- lapply(seq_along(ll), function(j) {
    res <- sapply(ll[[j]] , function(x) {
      A = x[[1]][coef.est != 0]$Gene
      B = x[[2]][coef.est != 0]$Gene
      length(intersect(A,B))/length(union(A,B))
    }) %>% mean
    names(res) <- paste0(clust_labs[[j]],"jacc")
    return(res)
  })

  clust_jacc %>% unlist
}

print("done clust and regress interaction")
