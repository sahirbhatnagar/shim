source("functions.R")
source("packages.R")

rm(list=ls())
p <- 1000
n <- 40
n0 <- 20
n1 <- n - n0
rho0 <- 0.1
rho1 <- 0.9



train <- simToyExample(p = 1000, n = 40, n0 = 20, rho0 = 0.1, rho1 = 0.9)
test <- simToyExample(p = 1000, n = 40, n0 = 20, rho0 = 0.1, rho1 = 0.9)

result <- generate_data_toy(p = 1000, Xtrain = train, Xtest = test,
                            mu_pc1 = c(rep(-2, 20), rep(2, 20)), 
                            mu_pc2 = c(rep(-1,10),rep(1,10),rep(-1,10),rep(1,10)),
                            alphaSim = 0.5, n = 40, n0 = 20,
                            cluster_distance = "corr",
                            include_interaction = F,
                            signal_to_noise_ratio = 1,
                            EclustAddDistance = "fisherScore",
                            clustMethod = "hclust",
                            cutMethod = "gap",
                            distanceMethod = "euclidean",
                            method = "ward.D2",
                            K.max = 10, B = 5)

# pheatmap::pheatmap(t(train), cluster_rows = F, cluster_cols = F, show_rownames = F)
# pheatmap::pheatmap(t(test), cluster_rows = F, cluster_cols = F, show_rownames = F)


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
                                    include_interaction = includeInteraction,
                                    s0 = result[["S0"]],
                                    p = p,
                                    gene_groups = result[["clustersAll"]],
                                    clust_type = "clust"),
                    SIMPLIFY = F,
                    USE.NAMES = F)