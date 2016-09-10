##################################
# R source code file for conducting simulations
# on hydra cluster
# Git: this is on the eclust repo, sim2-modules-mammouth branch
# Created by Sahir,  August 24, 2016
# Updated:  Sept 6, 2016
# Notes:
# Simulations being run from /mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/august2016simulation/linear
# This is the same simulation begin run on mammouth cluster
# The code needs to be changed because I am running mclapply on mammouth
# cluster, and need to run job arrays on hydra cluster
# to transfer files from local laptop when connected to LDI network at the Jewish:
# scp simulation.R sahir.bhatnagar@d1p-hydratm01:/mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/august2016simulation/linear
# qsub command from /mnt/GREENWOOD_SCRATCH/sahir.bhatnagar/august2016simulation/linear folder:
# for i in {1..24..1} ; do qsub -v index=$i sim-hydra.sh ; done
# to submit to a specific node:
# for i in {19..24..1} ; do qsub -v index=$i sim-hydra.sh -l nodes=d1p-hydraex04.ldi.lan:ppn=1; done
##################################

rm(list=ls())
source("packages.R")
source("functions.R")
options(digits = 4, scipen=999)

# source("/home/bhatnaga/coexpression/may2016simulation/sim2-modules-mammouth/packages.R")
# source("/home/bhatnaga/coexpression/may2016simulation/sim2-modules-mammouth/functions.R")

source(paste(Sys.getenv("PBS_O_WORKDIR"),"packages.R", sep="/"))
source(paste(Sys.getenv("PBS_O_WORKDIR"),"functions.R", sep="/"))

parametersDf <- expand.grid(rho = c(0.2,0.90),
                            p = c(5000),
                            SNR = c(0.2,1,2),
                            n = c(400), # this is the total train + test sample size
                            # nActive = c(300), # must be even because its being split among two modules
                            #n0 = 200,
                            cluster_distance = c("tom","corr"),
                            Ecluster_distance = c("difftom", "diffcorr"),
                            rhoOther = 0.6,
                            betaMean = c(1),
                            alphaMean = c(0.5,2),
                            betaE = 2,
                            includeInteraction = TRUE,
                            includeStability = TRUE,
                            distanceMethod = "euclidean",
                            clustMethod = "hclust",
                            #cutMethod = "gap",
                            cutMethod = "dynamic",
                            agglomerationMethod = "average",
                            K.max = 10, B = 10, stringsAsFactors = FALSE)

parametersDf <- transform(parametersDf, n0 = n/2, nActive = p*0.05)
parametersDf <- parametersDf[which(parametersDf$cluster_distance=="tom" & parametersDf$Ecluster_distance=="difftom" | 
                                     parametersDf$cluster_distance=="corr" & parametersDf$Ecluster_distance=="diffcorr"),]
parameterIndex <- as.numeric(as.character(commandArgs(trailingOnly = T)[1]))

# parameterIndex = 3
simulationParameters <- parametersDf[parameterIndex,, drop = F]

print(simulationParameters)

## ---- generate-data ----
message("generating data")

p <- simulationParameters[,"p"];
n <- simulationParameters[,"n"];
n0 <- simulationParameters[,"n0"];
SNR <- simulationParameters[,"SNR"]
n1 <- n - n0
cluster_distance <- simulationParameters[,"cluster_distance"]
Ecluster_distance <- simulationParameters[,"Ecluster_distance"]
rhoOther <- simulationParameters[,"rhoOther"];
betaMean <- simulationParameters[,"betaMean"];
betaE <- simulationParameters[,"betaE"]
alphaMean <- simulationParameters[,"alphaMean"];
rho <- simulationParameters[,"rho"];
nActive <- simulationParameters[,"nActive"];
includeInteraction <- simulationParameters[,"includeInteraction"]
includeStability <- simulationParameters[,"includeStability"]
distanceMethod <- simulationParameters[,"distanceMethod"]
clustMethod <- simulationParameters[,"clustMethod"]
cutMethod <- simulationParameters[,"cutMethod"]
agglomerationMethod <- simulationParameters[,"agglomerationMethod"]
K.max <- simulationParameters[,"K.max"]
B <- simulationParameters[,"B"]

# in this simulation its blocks 3 and 4 that are important
# leaveOut:  optional specification of modules that should be left out 
# of the simulation, that is their genes will be simulated as unrelated 
# ("grey"). This can be useful when simulating several sets, in some which a module 
# is present while in others it is absent.
d0 <- simModule(n = n0, p = p, rho = c(0,0), exposed = FALSE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.01,
                maxCor = 1,
                corPower = 1,
                propNegativeCor = 0.3,
                backgroundNoise = 0.5,
                signed = FALSE,
                leaveOut = 1:4)

d1 <- simModule(n = n1, p = p, rho = c(rho, rho), exposed = TRUE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.4,
                maxCor = 1,
                corPower = 0.3,
                propNegativeCor = 0.3,
                backgroundNoise = 0.5,
                signed = FALSE)

# these should be the same. if they arent, its because I removed the red and
# green modules from the E=0 group
truemodule0 <- d0$setLabels
t0 <- table(truemodule0)
truemodule1 <- d1$setLabels
t1 <- table(truemodule1)
table(truemodule0,truemodule1)

# Convert labels to colors for plotting
moduleColors <- labels2colors(truemodule1)
table(moduleColors, truemodule1)

X <- rbind(d0$datExpr, d1$datExpr) %>%
  magrittr::set_colnames(paste0("Gene", 1:p)) %>%
  magrittr::set_rownames(paste0("Subject",1:n))

dim(X)
# pheatmap(cor(X))
# pheatmap(cor(d1$datExpr))
# pheatmap(cor(d0$datExpr))
# pheatmap(cor(d1$datExpr)-cor(d0$datExpr))
# pheatmap(WGCNA::TOMsimilarityFromExpr(X))
# pheatmap(WGCNA::TOMsimilarityFromExpr(d1$datExpr))
# pheatmap(WGCNA::TOMsimilarityFromExpr(d1$datExpr)-WGCNA::TOMsimilarityFromExpr(d0$datExpr))

# betaMainEffect <- vector("double", length = p)
# betaMainInteractions <- vector("double", length = p)
# # first assign random uniform to every gene in cluster 3 and 4,
# # then randomly remove so that thers only nActive left
# betaMainEffect[which(truemodule1 %in% 3:4)] <- runif(sum(truemodule1 %in% 3:4),
#                                                      betaMean - 0.1, betaMean + 0.1)
#
# # randomly set some coefficients to 0 so that there are only nActive non zero
# betaMainEffect <- replace(betaMainEffect,
#                           sample(which(truemodule1 %in% 3:4), sum(truemodule1 %in% 3:4) - nActive,
#                                  replace = FALSE), 0)
#
# betaMainInteractions[which(betaMainEffect!=0)] <- runif(nActive, alphaMean - 0.1, alphaMean + 0.1)
#
# beta <- c(betaMainEffect,
#           betaE,
#           betaMainInteractions)
# plot(beta)

betaMainEffect <- vector("double", length = p)
betaMainInteractions <- vector("double", length = p)

# the first nActive/2 in the 3rd block are active
betaMainEffect[which(truemodule1 %in% 3)[1:(nActive/2)]] <- runif(
  nActive/2, betaMean - 0.1, betaMean + 0.1)

# the first nActive/2 in the 4th block are active
betaMainEffect[which(truemodule1 %in% 4)[1:(nActive/2)]] <- runif(
  nActive/2, betaMean - 0.1, betaMean + 0.1)

betaMainInteractions[which(betaMainEffect!=0)] <- runif(nActive, alphaMean - 0.1, alphaMean + 0.1)

# must be in this order!!!! main effects, E, and then interactions... this order is being used
# by the generate_data function
beta <- c(betaMainEffect,
          betaE,
          betaMainInteractions)

#plot(beta)

result <- generate_data(p = p, X = X, 
                        beta = beta, include_interaction = includeInteraction,
                        cluster_distance = cluster_distance,
                        n = n, n0 = n0, 
                        eclust_distance = Ecluster_distance,
                        signal_to_noise_ratio = SNR,
                        distance_method = distanceMethod,
                        cluster_method = clustMethod,
                        cut_method = cutMethod,
                        agglomeration_method = agglomerationMethod,
                        K.max = K.max, B = B, nPC = 1)

# result$clustersEclust[which(betaMainEffect!=0)][, table(module, cluster)]
# result$clustersEclust[module=="yellow"]
#
# result$clustersAll[which(betaMainEffect!=0)][, table(module, cluster)]
# result$clustersAll[module=="blue"]
## ---- univariate-pvalue -----

# message("Starting univariate p-value with interaction")
# 
# # filtering on univariate p-value, taking the lowet 5 percent and then running
# # multiple linear regression on it
# # the correlations are really messing up the oracle model which is why when
# # you filter you get better fit
# # stability = TRUE makes the function return coefficients only. use only for
# # calculating stability measures
# # output is in the following format: eg. uni_na_lm_yes_mse
# # which represents method_clustermeasure_model_includeE_measure
# 
# percent <- if (p > 1000) 0.005 else 0.01
# uni_res <- uniFit(train = result[["DT_train"]], test = result[["DT_test"]],
#                   percent = percent, stability = F,
#                   include_E = T,
#                   include_interaction = includeInteraction,
#                   filter_var = F, p = p,
#                   s0 = result[["S0"]], true_beta = result[["beta_truth"]])
# 
# if (includeStability) {
#   uni_stab <- mapply(uniFit,
#                      train = result[["DT_train_folds"]],
#                      MoreArgs = list(test = result[["DT_test"]],
#                                      s0 = result[["S0"]],
#                                      true_beta = result[["beta_truth"]],
#                                      stability = T,
#                                      include_E = T,
#                                      include_interaction = includeInteraction,
#                                      percent = percent,
#                                      filter_var = F,
#                                      p = p),
#                      SIMPLIFY = F)
#   
#   # Make the combinations of list elements
#   ll <- combn(uni_stab, 2, simplify = F)
#   
#   # Pairwise correlations of the model coefficients for each of the 10 CV folds
#   uni_mean_stab <- lapply(c("pearson","spearman"), function(i) {
#     res <- mean(sapply(ll,
#                        function(x) WGCNA::cor(x[[1]]$coef.est,
#                                               x[[2]]$coef.est,
#                                               method = i,use = 'pairwise.complete.obs')
#     ), na.rm = TRUE)
#     names(res) <- paste0("uni_na_lm_yes_",i)
#     return(res)
#   }
#   )
#   
#   # Jaccard index
#   uni_jacc <- mean(sapply(ll, function(x) {
#     A = x[[1]][coef.est != 0]$Gene
#     B = x[[2]][coef.est != 0]$Gene
#     if (length(A)==0 | length(B)==0) 0 else length(intersect(A,B))/length(union(A,B))
#   }
#   ), na.rm = TRUE)
#   
# }
# 
# message("done univariate p-value with interaction")
# 
# uni_res

## ---- cluster-and-regress ----

# we will treat the clusters as fixed i.e., even if we filter, or
# do cross validation, the group labels are predetermined by the
# above clustering procedure
# This method is based on clusters derived without accounting for the
# environment

print("starting cluster and regress with interaction")

clust_res <- mapply(clust_fun,
                    #summary = rep(c("pc","spc","avg"), each = 3),
                    #model = rep(c("lm", "lasso","elasticnet"), 3),
                    summary = rep(c("avg","pc"), each = 2),
                    model = rep(c("lasso","elasticnet"), 2),
                    MoreArgs = list(x_train = result[["X_train"]],
                                    x_test = result[["X_test"]],
                                    y_train = result[["Y_train"]],
                                    y_test = result[["Y_test"]],
                                    stability = F,
                                    filter = F,
                                    filter_var = F,
                                    include_E = T,
                                    include_interaction = includeInteraction,
                                    true_beta = result[["beta_truth"]],
                                    s0 = result[["S0"]],
                                    p = p,
                                    gene_groups = result[["clustersAll"]],
                                    clust_type = "clust",
                                    nPC = 1),
                    SIMPLIFY = F,
                    USE.NAMES = F)

# result %>% names
clust_res %>% unlist

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
                                                              include_interaction = includeInteraction,
                                                              gene_groups = result[["clustersAll"]],
                                                              p = p,
                                                              true_beta = result[["beta_truth"]],
                                                              clust_type = "clust",
                                                              nPC = 1),
                                              SIMPLIFY = F),
                       #summary = rep(c("pc","spc","avg"), each = 3),
                       #model = rep(c("lm", "lasso","elasticnet"), 3),
                       summary = rep(c("avg","pc"), each = 2),
                       model = rep(c("lasso","elasticnet"), 2),
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
                       summary = rep(c("avg","pc"), each = 2),
                       model = rep(c("lasso","elasticnet"), 2),
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

Eclust_res <- mapply(clust_fun,
                     #summary = rep(c("pc","spc","avg"), each = 3),
                     #model = rep(c("lm", "lasso","elasticnet"), 3),
                     summary = rep(c("avg","pc"), each = 2),
                     model = rep(c("lasso","elasticnet"), 2),
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
                                     true_beta = result[["beta_truth"]],
                                     gene_groups = result[["clustersAddon"]],
                                     clust_type = "Eclust",
                                     nPC = 1),
                     SIMPLIFY = F,
                     USE.NAMES = F)

# Eclust_res %>% names
# options(digits = 2, scipen=999)
# Eclust_res %>% unlist

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
                                                               true_beta = result[["beta_truth"]],
                                                               clust_type = "Eclust",
                                                               nPC = 1),
                                               SIMPLIFY = F),
                        #summary = rep(c("pc","spc","avg"), each = 3),
                        #model = rep(c("lm", "lasso","elasticnet"), 3),
                        summary = rep(c("avg","pc"), each = 2),
                        model = rep(c("lasso","elasticnet"), 2),
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
                        summary = rep(c("avg","pc"), each = 2),
                        model = rep(c("lasso","elasticnet"), 2),
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

## ---- group-lasso-clust ----

# currently not being used for simulation 2 because it takes way too long
# result[["gene_groups_inter"]] is without using Environment
# result[["gene_groups_inter_Addon"]] is with using Environment
# group_pen_res <- mapply(group_pen_fun, model = c("gglasso"),
#                         MoreArgs = list(x_train = result[["X_train"]], 
#                                         x_test = result[["X_test"]], 
#                                         y_train = result[["Y_train"]], 
#                                         y_test = result[["Y_test"]], 
#                                         stability = F,
#                                         filter = F,
#                                         filter_var = F,
#                                         include_E = T,
#                                         include_interaction = includeInteraction,
#                                         s0 = result[["S0"]],
#                                         p = p,
#                                         true_beta = result[["beta_truth"]],
#                                         gene_groups = result[["gene_groups_inter"]]),
#                         SIMPLIFY = F, 
#                         USE.NAMES = F)
# 
# group_pen_res %>% unlist()
# if (includeStability) {
#   group_pen_stab <- mapply(function(model) mapply(group_pen_fun, 
#                                                   x_train = result[["X_train_folds"]], 
#                                                   y_train = result[["Y_train_folds"]], 
#                                                   MoreArgs = list(stability = T, 
#                                                                   model = model,
#                                                                   filter = F,
#                                                                   filter_var = F,
#                                                                   include_E = T, 
#                                                                   include_interaction = T,
#                                                                   s0 = result[["S0"]],
#                                                                   true_beta = result[["beta_truth"]],
#                                                                   gene_groups = result[["gene_groups_inter"]],
#                                                                   p = 1000), 
#                                                   SIMPLIFY = F),
#                            model = c("gglasso"),
#                            SIMPLIFY = F, 
#                            USE.NAMES = F)
#   
#   
#   # Make the combinations of list elements
#   ll <- lapply(seq_along(group_pen_stab), function(i) combn(group_pen_stab[[i]], 2, simplify = F))
#   
#   
#   group_pen_labels <- function(model) {
#     paste0("group_pen",paste0("_",model),"_")
#   }
#   
#   group_pen_labs <- mapply(group_pen_labels, 
#                            model = c("gglasso"),
#                            USE.NAMES = F)
#   
#   
#   # Pairwise correlations of the model coefficients for each of the 10 CV folds 
#   group_pen_mean_stab <- lapply(seq_along(ll), function(j) {
#     lapply(c("pearson","spearman"), function(i) { 
#       res <- sapply(ll[[j]] , function(x) WGCNA::cor(x[[1]]$coef.est, x[[2]]$coef.est, method = i)) %>% mean
#       names(res) <- paste0(group_pen_labs[[j]], i)
#       return(res)
#     }
#     )
#   }
#   )
#   
#   group_pen_mean_stab %>% unlist
#   
#   # Jaccard index
#   group_pen_jacc <- lapply(seq_along(ll), function(j) {
#     res <- sapply(ll[[j]] , function(x) {  
#       A = x[[1]][coef.est != 0]$Gene
#       B = x[[2]][coef.est != 0]$Gene
#       length(intersect(A,B))/length(union(A,B))
#     }) %>% mean
#     names(res) <- paste0(group_pen_labs[[j]],"jacc")
#     return(res)
#   })
#   
#   group_pen_jacc %>% unlist
# }
# print("done group penalization no interaction")



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

## ---- final-results ----

final_results <- if (includeStability) {
  c(simulationParameters,
    # uni_res,uni_mean_stab, "uni_jacc" = uni_jacc,
    clust_res,clust_mean_stab, clust_jacc,
    pen_res,pen_mean_stab, pen_jacc,
    Eclust_res,Eclust_mean_stab, Eclust_jacc) %>% unlist } else {
      c(simulationParameters, 
        # uni_res,
        clust_res,
        pen_res,
        Eclust_res)  %>% unlist
    }

final_results %>% t %>% as.data.frame()


filename <- tempfile(pattern = paste0(sprintf("%s_%s_rho%.2f_p%1.0f_SNR%.2f_n%1.0f_s0%1.0f_beta%.2f_alpha%.2f",
                                              cluster_distance,Ecluster_distance,rho,p,SNR, n, nActive, betaMean, alphaMean),"_"),
                     tmpdir = paste(Sys.getenv("PBS_O_WORKDIR"),"results", sep="/"))
#tmpdir = "/home/bhatnaga/coexpression/may2016simulation/sim2-modules-mammouth/results/")
write.table(final_results %>% t %>% as.data.frame(),
            file = filename,
            quote = F,
            row.names = F,
            col.names = F)

write.table(final_results %>% t %>% as.data.frame() %>% colnames(),
            file = paste(Sys.getenv("PBS_O_WORKDIR"),"colnames_stab_hydra.txt",sep="/"),
            # file  = filename,
            quote = F,
            row.names = F, col.names = F)

# final_results %>% t %>% as.data.frame() %>% colnames()


