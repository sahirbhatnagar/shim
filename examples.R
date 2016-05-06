# # number of observations
# n <- 100
#
# # number of predictors
# p <- 5
#
# # environment variable
# e <- sample(c(0,1), n, replace = T)
#
# # main effects
# x <- cbind(matrix(rnorm(n*p), ncol = p), e)
#
# # need to label columns
# dimnames(x)[[2]] <- c("cg12456","cg45","cg88","cg7888","cg87979","environ")
#
# # design matrix without intercept
# X <- model.matrix(~e*(x1+x2+x3+x4+x5)-1, data = as.data.frame(x))
# X <- model.matrix(~environ*(cg12456+cg45+cg88)+cg7888*cg45+cg87979*cg7888-1, data = as.data.frame(x))
# head(X)
#
# interaction_names <- grep(":", colnames(X), value = T)
# main_effect_names <- setdiff(colnames(X), interaction_names)
#
#
# # response
# Y <- X %*% rbinom(ncol(X), 1, 0.6) + 3*rnorm(n)
#
# # standardize data
# data_std <- standardize(X,Y)
# head(data_std$x);head(data_std$y)
# ridge_weights(x = data_std$x, y = data_std$y,
#               main.effect.names = main_effect_names,
#               interaction.names = interaction_names)
#
# shim_once(x = data_std$x, y = data_std$y,
#           main.effect.names = main_effect_names,
#           interaction.names = interaction_names)
#
# rm(list = ls())


# example 1 ---------------------------------------------------------------
library(magrittr)
library(data.table)
library(dplyr)
library(eclust)
rm(list=ls())
set.seed(123456)
# number of observations
n <- 1000

# number of predictors
p <- 5

# environment variable
e <- sample(c(0,1), n, replace = T)

# main effects
x <- cbind(matrix(rnorm(n*p), ncol = p), e)

# need to label columns
dimnames(x)[[2]] <- c("x1","x2","x3","x4","x5","e")

# design matrix without intercept (can be user defined interactions)
X <- model.matrix(~(x1+x2+x3)*e+x1*x4+x3*x5-1, data = as.data.frame(x))
head(X)
# names must appear in the same order as X matrix
interaction_names <- grep(":", colnames(X), value = T)
main_effect_names <- setdiff(colnames(X), interaction_names)

# response
(trueBeta <- c(0,0,2,1.5,0,0,0,0,4,0,0))
Y <- X %*% trueBeta + 3*rnorm(n)


res <- shim(x = X, y = Y,
     main.effect.names = main_effect_names,
     interaction.names = interaction_names,
     initialization.type = "univariate")

library(doMC)
registerDoMC(cores = 3)
cvres <- cv.shim(x = X, y = Y,
            main.effect.names = main_effect_names,
            interaction.names = interaction_names,
            initialization.type = "ridge",
            nlambda.gamma = 5, nlambda.beta = 5, nlambda = 25,
            max.iter = 50, threshold = 1e-3)

names(cvres)
class(cvres)
library(ggplot2)
library(latex2exp)
plot(cvres)
coef(cvres, s="lambda.min")
coef(cvres, s="lambda.1se")
names(cvres)
cvres$lambda.1se.beta
cvres$shim.fit$call
class(cvres$shim.fit)
coef(cvres$shim.fit)
names(cvres$shim.fit)

res$tuning.parameters[,c("s3","s4"),drop=F]
res$beta[,c("s3","s4"),drop=F]
res$alpha[,c("s3","s4"),drop=F]

plot(log(res$tuning.parameters["lambda.gamma",]),
          log(res$lambda.gamma))

options(digits = 4, scipen = 999)
resSingle <- shim(x = X, y = Y,
            main.effect.names = main_effect_names,
            interaction.names = interaction_names,
            lambda.beta = res$lambda.beta[1:10],
            lambda.gamma = res$lambda.gamma[1:10],
            nlambda = 10, nlambda.gamma = 100, nlambda.beta = 100,
            verbose = T, initialization.type = "univariate")

res$tuning.parameters[,1:10]
res$beta[,1:5]
resSingle$beta[,1:5]
res$alpha[,1:5]
resSingle$alpha[,1:5]

res$beta[,1:15]
resSingle$beta
res$gamma[,"s"]

resSingle$beta[,"s3"]
resSingle$gamma[,"s3"]

resSingle$beta
resSingle$alpha

resSingle$beta
resSingle$alpha




args(shim)

(nd <- dim(res$tuning.parameters[, c("s1","s2")]))
(nd <- dim(res$tuning.parameters[, c("s1")]))
(nd <- dim(res$tuning.parameters["lambda.beta","s2"]))
if(FALSE | nd[2]>1) 1
res$tuning.parameters["lambda.beta","s2", drop=F] %>% length
res$tuning.parameters["lambda.gamma","s2", drop=F] %>% length

res$call
shim

knitr::kable(as.matrix(res$beta))
DT::datatable(round(res$tuning.parameters,5))
names(res)

res$tuning.parameters %>% str
res$dfbeta
res$dfalpha
eclust::nonzero(res$beta)

require(doMC)
registerDoMC(cores = 10)
system.time(cv.res <- cv.shim(x = data_std$x, y = data_std$y,
            main.effect.names = main_effect_names,
            interaction.names = interaction_names))

require(ggplot2)
require(latex2exp)
cv.res %>% names
plot(cv.res)

cv.res$lambda.min.name
cv.res$lambda.1se.name
cv.res$lambda.

# example 2 ---------------------------------------------------------------
devtools::load_all()
library(magrittr)
library(data.table)
library(dplyr)
rm(list=ls())
set.seed(123456)

# number of predictors
p = 10

# number of test subjects
n = 200

# correlation between X's
rho = 0.45

# signal to noise ratio
signal_to_noise_ratio = 4

# names of the main effects, this will be used in many of the functions
main_effect_names <- paste0("x",1:p)

# names of the active set
true_var_names <- c("x1","x2","x3","x4","x1:x2", "x1:x3", "x1:x4", "x2:x3", "x2:x4", "x3:x4")

# different true coefficient vectors as in Table 1 of Choi et al.
beta1 <- c(7,2,1,1,0,0,0,0,0,0) %>% magrittr::set_names(true_var_names)
beta2 <- c(7,2,1,1,1,0,0,0.5,0.4,0.1) %>% magrittr::set_names(true_var_names)
beta3 <- c(7,2,1,1,7,7,7,2,2,1) %>% magrittr::set_names(true_var_names)
beta4 <- c(7,2,1,1,14,14,14,4,4,2) %>% magrittr::set_names(true_var_names)
beta5 <- c(0,0,0,0,7,7,7,2,2,1) %>% magrittr::set_names(true_var_names)

# simulate Toeplitz like correlation structure between X's
H <- abs(outer(1:p, 1:p, "-"))
cor <- rho^H

# generate X's from multivariate normal and label the matrix
DT <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = cor) %>%
  magrittr::set_colnames(paste0("x",1:p)) %>%
  set_rownames(paste0("Subject",1:n))

# create X matrix which contains all main effects and interactions
# but not the intercept
# each column is standardized to mean 0 and sd 1
X <- model.matrix(
  as.formula(paste0("~(",paste0(main_effect_names, collapse = "+"),")^2-1")),
  data = DT %>% as.data.frame()) # %>% scale

# not doing this before.. now handled by shim_multiple_faster function
# check that means of columns are 0 and sd 1
#colMeans(X) %>%  sum
#apply(X, 2, sd) %>% sum

# generate response with user defined signal to noise ratio and center
# the response
y.star <- X[,names(beta3)] %*% beta3
error <- rnorm(n)
k <- sqrt(var(y.star)/(signal_to_noise_ratio*var(error)))
Y <- y.star + k*error # %>% scale(center = TRUE, scale = FALSE)
colnames(Y) <- "Y"

# record mean of response before centering
(b0 <- mean(y.star + k*error))

# names of interaction variables assuming interaction terms contain a ":"
# this will be used in many of the functions
interaction_names <- colnames(X) %>% grep(":",., value = T)


# example 3 ---------------------------------------------------------------

library(magrittr)
library(data.table)
library(dplyr)
rm(list=ls())
source("https://raw.githubusercontent.com/noamross/noamtools/master/R/proftable.R")
set.seed(123456)
# number of observations
n <- 100

# number of predictors
p <- 10

# correlation between X's
rho <- 0.75

# signal to noise ratio
signal_to_noise_ratio <- 4

# environment variable
e <- sample(c(0,1), n, replace = T)

# simulate Toeplitz like correlation structure between X's
H <- abs(outer(1:p, 1:p, "-"))

# # constant correlation structure
# H <- abs(outer(rep(2,p), rep(1,p), "-"))
# diag(H) <- 0
cor <- rho^H

# generate X's from multivariate normal and label the matrix
DT <- cbind(MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = cor), e) %>%
  magrittr::set_colnames(c(paste0("x",seq_len(p)),"e")) %>%
  set_rownames(paste0("Subject",1:n))

head(DT)

# design matrix without intercept (can be user defined interactions)
(form <- paste0("~(",paste(paste0("x",seq_len(p)), collapse = "+"),")*e - 1" ))

X <- model.matrix(as.formula(form), data = as.data.frame(DT))
head(X)

# response
# names of the active set
(true_var_names <- c(paste0("x",1:5),"e",paste0("x",1:5,":e") ))
trueBeta <- rep(0, ncol(X))
trueBeta[which(colnames(X) %in% true_var_names)] <- c(runif(5,3.9,4.1),4, runif(5,1.9,2.1))
trueBeta

# generate response with user defined signal to noise ratio and center
# the response
y.star <- X %*% trueBeta
error <- rnorm(n)
k <- sqrt(var(y.star)/(signal_to_noise_ratio*var(error)))
Y <- y.star + k*error
colnames(Y) <- "Y"

# record mean of response before centering
(b0 <- mean(y.star + k*error))

# names of interaction variables assuming interaction terms contain a ":"
# this will be used in many of the functions
# names must appear in the same order as X matrix
(interaction_names <- grep(":", colnames(X), value = T))
(main_effect_names <- setdiff(colnames(X), interaction_names))

coef(lm.fit(x = x[, 1, drop = F], y = y)) %>% magrittr::extract(1)
coef(lm.fit(x = x[, 1, drop = F], y = y))[1]
summary(lm(Y~., data = as.data.frame(cbind(Y,X))))


res <- shim(x = X, y = Y,
            main.effect.names = main_effect_names,
            interaction.names = interaction_names)
require(doMC)
registerDoMC(cores = 3)
system.time(cv.res <- cv.shim(x = X, y = Y,
                              main.effect.names = main_effect_names,
                              interaction.names = interaction_names,
                              parallel = TRUE,
                              nfolds = 5,
                              nlambda.beta = 100,
                              nlambda.gamma = 1,
                              nlambda = 100,
                              verbose  = FALSE))
library(latex2exp)
library(ggplot2)
plot(cv.res)
coef(cv.res, s = "lambda.min")


names(cv.res)
cv.res$shim.fit
