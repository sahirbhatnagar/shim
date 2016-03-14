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


library(magrittr)
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

library(dplyr)
coef(cv.res)
coef.cv.shim
shim_once(x = data_std$x, y = data_std$y,
          main.effect.names = main_effect_names,
          interaction.names = interaction_names,
          nlambda.gamma = 5, nlambda.beta = 5)


DT <- data.frame(y = 1:nobs, x = ftime[ord], group = sample(c("E=1","E=0"),nobs, replace = T))
str(DT)

X <- data.frame(
  y = c(1,1,nobs,nobs),
  x = c(0,0,0,0),
  group = c("E=1","E=0","E=1","E=0"))

dat <- rbind(DT,X)

dat[]
order()
p <- ggplot(dat, aes(x=x, y=y))
p + geom_polygon() + facet_grid(~group)





chull(X)
## Not run:
# Example usage from graphics package
plot(X, cex = 0.5)
hpts <- chull(X)
hpts <- c(hpts, hpts[1])
lines(X[hpts, ])



library(survival)
data(veteran)
table(veteran$status)
evtimes <- veteran$time[veteran$status == 1]
hist(evtimes, nclass=30, main='', xlab='Survival time (days)', col='gray90', probability=TRUE)
tgrid <- seq(0, 1000, by=10)
lines(tgrid, dexp(tgrid, rate=1.0/mean(evtimes)),
      lwd=2, lty=2, col='red')

veteran$prior <- factor(veteran$prior, levels = c(0, 10))
veteran$celltype <- factor(veteran$celltype,
                           levels = c('large', 'squamous', 'smallcell', 'adeno'))
veteran$trt <- factor(veteran$trt, levels = c(1, 2))

y <- with(veteran, Surv(time, status))
nobs <- nrow(y)
ftime <- veteran$time
ord <- order(ftime, decreasing=TRUE)
plot(0, type='n', xlim=c(0, max(ftime)), ylim=c(0, nobs),
     xlab='Follow-up time', ylab='Population')
segments(rep(0.0, nobs), 1:nobs, ftime[ord], 1:nobs, col='gray25')
cases <- veteran$status == 1
points((ftime[ord])[cases[ord]], (1:nobs)[cases[ord]], pch=20, col='red', cex=0.5)



# Simulate censored survival data for two outcome types from Weibull distributions:

set.seed(1)
nobs <- 5000

a1 <- 1.0
b1 <- 200
a2 <- 1.0
b2 <- 50
c1 <- 0.0
c2 <- 0.0

z <- rbinom(nobs, 1, 0.5)
t1 <- rweibull(nobs, a1, b1 * exp(z * c1)^(-1/a1))
t2 <- rweibull(nobs, a2, b2 * exp(z * c2)^(-1/a2))
tlim <- 10
mean(t1 < tlim)
mean(t2 < tlim)
e <- 1 * (t1 < t2) + 2 * (t1 >= t2)
t <- pmin(t1, t2)
e[t >= tlim] <- 0
t[t >= tlim] <- tlim
table(e)

evtimes=t[e==1]
ncc <- length(unique(e))

tgrid=sort(unique(c(seq(0, max(t), 0.1), evtimes)))
nt <- length(tgrid)
ne <- matrix(NA, ncc, nt)
dim(nriskset)
for (i in 1:nt) {
  for (j in 1:ncc) {
    ne[j,i]= sum(t[e==(j-1)] < tgrid[i])
  }
}

y = ne[2,match(evtimes, tgrid)] + (nobs - (colSums(ne[,match(evtimes, tgrid)]))) * runif(sum(e==1))

N <- nobs
max.yr=max(t)
dev.off()
plot(1,1, xlab="Time",ylab="Population",ylim=c(0,N),xlim=c(0,max.yr),type="n",cex.lab=1.2)
mtext("Population", side = 2, cex=1.2,line=4)
rect(0,0,max.yr,nobs,col="grey90",border=NA)

# Competing events remove individuals from the riskset, presented in the cyan area:
polygon(c(tgrid,max.yr,max.yr,0),c(nobs-(ne[1,]+ne[3,]),nobs-max(ne[1,]+ne[3,]),nobs,nobs),col='cyan',border=NA)

# If random censoring, this can be presented with a different colour:
#polygon(c(tgrid,max.yr,max.yr,0),c(nobs-ne[1,],nobs-max(ne[1,]),nobs,nobs),col='black',border=NA)

# Individuals removed from the risksets due to events of interest presented in the red area:
polygon(c(tgrid,max.yr,max.yr,0),c(ne[2,],max(ne[2,]),0,0),col='red',border=NA)

# Events of interest:
points(evtimes, y, col='red', pch=16, cex=.5)




