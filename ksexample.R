source("~/git_repositories/eclust/k_filter.R")

n<-200
p<-5000

x<-matrix(runif(n*p),n,p)
#x<-x%*%sigma.sqrt

y <- 4*x[,1]+2*tan(x[,2]*pi/2)+5*x[,3]^2+rnorm(n)
obj <- k.filter(x=x,y=y)

# the minimum number of predictors needed to keep all useful predictors
# in this case, all useful predictor are x1,x2 and x3
max(obj$k.rank[c(1:3)])



obj$k.rank



# each row probably represents the number of slices used to categorize the
# response. each column of k.stat.single is a predictor
dim(obj$k.stat.single)
obj$k.stat.single[,1:10]


obj$k.stat



obj2 <- k.filter(x=X,y=Y)

obj2$k.rank
