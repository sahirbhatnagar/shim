# number of observations
n <- 100

# number of predictors
p <- 5

# environment variable
e <- sample(c(0,1), n, replace = T)

# main effects
x <- cbind(matrix(rnorm(n*p), ncol = p), e)

# need to label columns
dimnames(x)[[2]] <- c("cg12456","cg45","cg88","cg7888","cg87979","environ")

# design matrix without intercept
X <- model.matrix(~e*(x1+x2+x3+x4+x5)-1, data = as.data.frame(x))
X <- model.matrix(~environ*(cg12456+cg45+cg88)+cg7888*cg45+cg87979*cg7888-1, data = as.data.frame(x))
head(X)

interaction_names <- grep(":", colnames(X), value = T)
main_effect_names <- setdiff(colnames(X), interaction_names)


# response
Y <- X %*% rbinom(ncol(X), 1, 0.6) + 3*rnorm(n)

# standardize data
data_std <- standardize(X,Y)
head(data_std$x);head(data_std$y)
ridge_weights(x = data_std$x, y = data_std$y,
              main.effect.names = main_effect_names,
              interaction.names = interaction_names)

shim_once(x = data_std$x, y = data_std$y,
          main.effect.names = main_effect_names,
          interaction.names = interaction_names)



# number of observations
n <- 100

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

# names must appear in the same order as X matrix
interaction_names <- grep(":", colnames(X), value = T)
main_effect_names <- setdiff(colnames(X), interaction_names)

# response
Y <- X %*% rbinom(ncol(X), 1, 0.6) + 3*rnorm(n)

# standardize data
data_std <- standardize(X,Y)


res <- shim(x = data_std$x, y = data_std$y,
     main.effect.names = main_effect_names,
     interaction.names = interaction_names)


res$alpha


res$tuning.parameters %>% str
res$dfbeta
res$dfalpha
eclust::nonzero(res$beta)

require(doMC)
registerDoMC(cores = 10)
system.time(cv.res <- cv.shim(x = data_std$x, y = data_std$y,
            main.effect.names = main_effect_names,
            interaction.names = interaction_names))

cv.res %>% names


library(dplyr)
coef(cv.res)
coef.cv.shim
shim_once(x = data_std$x, y = data_std$y,
          main.effect.names = main_effect_names,
          interaction.names = interaction_names,
          nlambda.gamma = 5, nlambda.beta = 5)



