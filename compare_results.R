##################################
# R source code file for analyzing bouchard data and comparing Rsquared between different methods
# Created by Sahir, September 26, 2016
# Updated:
# NOTE: Hosted on 'bouchard' branch of eclust repo on github.
# see also 'clust_kmeans_bouchard_qsub.R'
# 
##################################



# Load data ---------------------------------------------------------------

rm(list=ls())
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/functions.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/packages.R")
source("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/data_cleaning.R")
load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/1_PC_bouchard_sd_filter_5K_probes_TOM_DIFFTOM_Sept19.RData")
res5k <- pp
rm(pp)
load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/1_PC_bouchard_sd_filter_10K_probes_TOM_DIFFTOM_Sept19.RData")
res10k <- pp
rm(pp)
load("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/bouchard_git/1_PC_bouchard_sd_filter_20K_probes_TOM_DIFFTOM_Sept19.RData")
res20k <- pp
rm(pp)

probe_sd <- rowSds(placentaALL)

# 5,000 most variable probes
filterd_probes <- probe_sd[probe_sd > quantile(probe_sd,0.978245)] %>% names
filterd_probes %>% length()

# 10,000 most variable probes
# filterd_probes <- probe_sd[probe_sd > quantile(probe_sd,0.956485)] %>% names
# filterd_probes %>% length()

# Create training and test set indices ------------------------------------

Y <- DT.pheno.placenta$`imc z-score` %>% as.matrix() %>% drop

trainIndex <- caret::createDataPartition(Y, p = 1, list = FALSE, times = 1) %>% drop

# not really used... need to do LOOCV 
# testIndex <- which(seq_len(length(Y)) %ni% trainIndex)
testIndex <- trainIndex
# all(c(trainIndex, testIndex) %in% seq_len(length(Y)))

# we also need the exposure variable as a vector
# columns of methylation data match the ordering of the phenotype data
# see 'data_cleaning.R' for details
# needs to be 0 and 1
E <- as.numeric(DT.pheno.placenta$case)-1


# Different input data (all based on TOM, not using correlations, no test set) ----

# 5k most variable 
pc_train_5k_eclust <- res5k$clustersAddon$PC
pc_train_5k_clust <- res5k$clustersAll$PC
avg_train_5k_eclust <- res5k$clustersAddon$averageExpr
avg_train_5k_clust <- res5k$clustersAll$averageExpr
original_train_5k <- t(placentaALL[filterd_probes,])

# eclust and clust interaction data 
pc_train_5k_eclust_interaction <- prepare_data(data = cbind(pc_train_5k_eclust, 
                                                            pheno = Y[trainIndex],
                                                            income = E[trainIndex]),
                                               response = "pheno", exposure = "income")

pc_train_5k_clust_interaction <- prepare_data(data = cbind(pc_train_5k_clust, 
                                                            pheno = Y[trainIndex],
                                                            income = E[trainIndex]),
                                               response = "pheno", exposure = "income")

avg_train_5k_eclust_interaction <- prepare_data(data = cbind(avg_train_5k_eclust, 
                                                            pheno = Y[trainIndex],
                                                            income = E[trainIndex]),
                                               response = "pheno", exposure = "income")

avg_train_5k_clust_interaction <- prepare_data(data = cbind(avg_train_5k_clust, 
                                                            pheno = Y[trainIndex],
                                                            income = E[trainIndex]),
                                               response = "pheno", exposure = "income")

# original interaction data
original_train_5k_interaction <- prepare_data(data = cbind(as.data.frame(original_train_5k), 
                                                            pheno = Y[trainIndex],
                                                            income = E[trainIndex]),
                                               response = "pheno", exposure = "income")



# Fit models --------------------------------------------------------------

fitControl <-  trainControl(method = "repeatedcv",
                            number = 5,
                            repeats = 3,
                            verboseIter = TRUE)

fitControl <-  trainControl(method = "boot",
                            number = 25,
                            # repeats = 3,
                            verboseIter = TRUE)

# Define the candidate models to test
marsGrid <- expand.grid(.degree = 1:3, .nprune = 1000)

#############################################################################
#               MARS Original                                            ----
#############################################################################
set.seed(1056)
mars_original_5k <- train(original_train_5k_interaction$X[,original_train_5k_interaction$main_effect_names], 
                    original_train_5k_interaction$Y,
                    method = "earth",
                    trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                    tuneGrid = marsGrid,
                    trControl = fitControl)

#############################################################################
#               MARS Clust                                              ----
#############################################################################
set.seed(1056)
# marsGridClust <- expand.grid(.degree = 1, .nprune = 1000)
mars_pc_5k_clust <- train(pc_train_5k_clust_interaction$X[, pc_train_5k_clust_interaction$main_effect_names], 
                          pc_train_5k_clust_interaction$Y,
                          method = "earth",
                          trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                          tuneGrid = marsGrid,
                          trControl = fitControl)

set.seed(1056)
mars_avg_5k_clust <- train(avg_train_5k_clust_interaction$X[, avg_train_5k_clust_interaction$main_effect_names], 
                           avg_train_5k_clust_interaction$Y,
                           method = "earth",
                           trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                           tuneGrid = marsGrid,
                           trControl = fitControl)

#############################################################################
#                     MARS EClust                                        ----
#############################################################################
set.seed(1056)
mars_pc_5k_eclust <- train(pc_train_5k_eclust_interaction$X[, pc_train_5k_eclust_interaction$main_effect_names], 
                           pc_train_5k_eclust_interaction$Y,
                           method = "earth",
                           trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                           # Explicitly declare the candidate models to test
                           tuneGrid = marsGrid,
                           trControl = fitControl)

set.seed(1056)
mars_avg_5k_eclust <- train(avg_train_5k_eclust_interaction$X[, avg_train_5k_eclust_interaction$main_effect_names], 
                           avg_train_5k_eclust_interaction$Y,
                           method = "earth",
                           trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                           # Explicitly declare the candidate models to test
                           tuneGrid = marsGrid,
                           trControl = fitControl)


#############################################################################
#               Bagged MARS Original                                    ----
#############################################################################
set.seed(1056)
bag_mars_original_5k <- train(original_train_5k_interaction$X[,original_train_5k_interaction$main_effect_names], 
                              original_train_5k_interaction$Y,
                              method = "bagEarth", B = 50,
                              trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                              tuneGrid = marsGrid,
                              trControl = fitControl)


#############################################################################
#               Bagged MARS Clust                                        ----
#############################################################################
set.seed(1056)
bag_mars_pc_5k_clust <- train(pc_train_5k_clust_interaction$X[, pc_train_5k_clust_interaction$main_effect_names], 
                          pc_train_5k_clust_interaction$Y,
                          method = "bagEarth", B = 50,
                          trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                          # Explicitly declare the candidate models to test
                          tuneGrid = marsGrid,
                          trControl = fitControl)

set.seed(1056)
bag_mars_avg_5k_clust <- train(avg_train_5k_clust_interaction$X[, avg_train_5k_clust_interaction$main_effect_names], 
                           avg_train_5k_clust_interaction$Y,
                           method = "bagEarth", B = 50,
                           trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                           # Explicitly declare the candidate models to test
                           tuneGrid = marsGrid,
                           trControl = fitControl)

#############################################################################
#               Bagged MARS EClust                                        ----
#############################################################################
set.seed(1056)
bag_mars_pc_5k_eclust <- train(pc_train_5k_eclust_interaction$X[, pc_train_5k_eclust_interaction$main_effect_names], 
                           pc_train_5k_eclust_interaction$Y,
                           method = "bagEarth", B = 50,
                           trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                           # Explicitly declare the candidate models to test
                           tuneGrid = marsGrid,
                           trControl = fitControl)

set.seed(1056)
bag_mars_avg_5k_eclust <- train(avg_train_5k_eclust_interaction$X[, avg_train_5k_eclust_interaction$main_effect_names], 
                            avg_train_5k_eclust_interaction$Y,
                            method = "bagEarth", B = 50,
                            trace = 4, nk = 1000, keepxy = TRUE, pmethod = "backward",
                            # Explicitly declare the candidate models to test
                            tuneGrid = marsGrid,
                            trControl = fitControl)


#############################################################################
#               Lasso Original                                            ----
#############################################################################

# Lasso + ECLUST ----
# apparently this does tune over both alpha and lambda, but just need to provide lambda
# http://stats.stackexchange.com/questions/69638/does-caret-train-function-for-glmnet-cross-validate-for-both-alpha-and-lambda?rq=1
# lassoGrid <- expand.grid(.alpha = seq(0,1,length.out = 5), .lambda = seq(min(glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda),
#                                                               max(glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda),length.out = 5))
# glmnet(x = kl_5k$X, y = kl_5k$Y)$lambda

set.seed(1056)
lasso_original_5k <- train(x = original_train_5k_interaction$X, 
                           y = original_train_5k_interaction$Y,
                           method = "glmnet",
                           tuneLength = 20,
                           trControl = fitControl)

#############################################################################
#               Lasso Clust                                            ----
#############################################################################

set.seed(1056)
lasso_pc_5k_clust <- train(x = pc_train_5k_clust_interaction$X, 
                            y = pc_train_5k_clust_interaction$Y,
                            method = "glmnet",
                            tuneLength = 20,
                            trControl = fitControl)

set.seed(1056)
lasso_avg_5k_clust <- train(x = avg_train_5k_clust_interaction$X, 
                            y = avg_train_5k_clust_interaction$Y,
                            method = "glmnet",
                            tuneLength = 20,
                            trControl = fitControl)

#############################################################################
#               Lasso EClust                                            ----
#############################################################################

set.seed(1056)
lasso_pc_5k_eclust <- train(x = pc_train_5k_eclust_interaction$X, 
                            y = pc_train_5k_eclust_interaction$Y,
                            method = "glmnet",
                            tuneLength = 20,
                            trControl = fitControl)

set.seed(1056)
lasso_avg_5k_eclust <- train(x = avg_train_5k_eclust_interaction$X, 
                             y = avg_train_5k_eclust_interaction$Y,
                             method = "glmnet",
                             tuneLength = 20,
                             trControl = fitControl)




resamps <- resamples(list(mars_original_5k=mars_original_5k,
                          mars_pc_5k_clust=mars_pc_5k_clust,
                          mars_avg_5k_clust=mars_avg_5k_clust,
                          mars_pc_5k_eclust=mars_pc_5k_eclust,
                          mars_avg_5k_eclust=mars_avg_5k_eclust,
                          bag_mars_pc_5k_clust=bag_mars_pc_5k_clust,
                          bag_mars_avg_5k_clust=bag_mars_avg_5k_clust,
                          bag_mars_pc_5k_eclust=bag_mars_pc_5k_eclust,
                          bag_mars_avg_5k_eclust=bag_mars_avg_5k_eclust))
                          # lasso_original_5k=lasso_original_5k,
                          # lasso_pc_5k_clust=lasso_pc_5k_clust,
                          # lasso_avg_5k_clust=lasso_avg_5k_clust,
                          # lasso_pc_5k_eclust=lasso_pc_5k_eclust,
                          # lasso_avg_5k_eclust=lasso_avg_5k_eclust))

resamps %>% str

dev.off()
trellis.par.set("theme1")
bwplot(resamps, horizontal=T, as.table = T)

dotplot(resamps, scales = list(x = list(relation = "free")), xlabel = "")

caret::densityplot.rfe(resamps)


summary(resamps)

caret:::bwplot.resamples(resamps, scales =list(x = list(relation = "free")))
splom(resamps)
densityplot(resamps)
xyplot(resamps)
caret:::densityplot.resamples(resamps, auto.key = list(columns = 3),
                              pch = "|")








