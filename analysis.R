##################################
# R source code file for analyzing NIHPD data
# Created by Sahir, May 24
# Updated:
# NOTE: This directory is more upto data than bouchard/scripts/ccareport
# I am now using version control in bouchard/cca/
# because it was getting too messy in ccareport
##################################

rm(list=ls())
source("functions.R")
source("packages.R")
source("data_cleaning.R")

brain_probes[1:10]
colnames(DT_with_pheno)[1:40]

DT_with_pheno[,"Subject_Gender",with=F] %>% str
DT_with_pheno[is.na(income)] %>% dim
DT_with_pheno[is.na(WASI_Full_Scale_IQ)] %>% dim
DT_with_pheno[is.na(Subject_Gender)] %>% dim
DT_with_pheno[is.na(Site_Location)] %>% dim


## ---- test-stat ----
registerDoMC(cores = 30)
res <- mclapply(brain_probes, function(i) {
  
  # i = brain_probes[1]
  dat <- DT_with_pheno[,c("WASI_Full_Scale_IQ","Subject_Gender","Site_Location",i),with=F]
  dat <- dat[!is.na(WASI_Full_Scale_IQ)]
  include_E <- FALSE
  include_interaction <- FALSE
  fit <- lm(as.formula(paste("WASI_Full_Scale_IQ ~",i, if (include_E) "+E", if (include_interaction) paste0("+",i,":E"),
                             "+ Subject_Gender + Site_Location")),
            data = dat)
  #fit %>% summary
  coef.index <- if (include_interaction) paste0(i,":E") else i
  
  data.frame(probe = i,
             pvalue = as.numeric(pt(abs(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
                                    df = fit$df.residual, lower.tail = F)*2),
             test.stat = as.numeric(coef(fit)[coef.index]/vcov(fit)[coef.index,coef.index] ^ 0.5),
             'mean' = mean(dat[,i, with = F][[1]]),
             'sd' = sd(dat[,i,with = F][[1]]), stringsAsFactors = FALSE)
}, mc.cores = 30)

uni_results <- do.call(rbind, res)
uni_results <- as.data.table(uni_results)


# uni_results[order(pvalue)[1:1000]]
# uni_results[, hist(sd)]
# uni_results[, hist(test.stat)]



## ---- filter-stats ----

with(uni_results,
     plot(rank(mean)/length(mean), -log10(pvalue), pch=16, cex=0.45))

par(mfrow = c(1,2))
with(uni_results,
     plot(rank(sd)/length(sd), -log10(pvalue), pch=16, cex=0.45))

trsf = function(n) log10(n+1)
plot(ecdf(trsf(uni_results$sd)), xlab=body(trsf), main="")
dev.off()


theta = seq(from=0, to=0.5, by=0.1)
pBH = filtered_p(filter=uni_results$mean, test=uni_results$pvalue, theta=theta, method="BH")
head(pBH)

rejection_plot(pBH, at="sample",
               xlim=c(0, 0.5), ylim=c(0, 15000),
               xlab="FDR cutoff (Benjamini & Hochberg adjusted p-value)", main="")

theta = seq(from=0, to=0.8, by=0.02)
rejBH = filtered_R(alpha = 0.1,
                   filter=uni_results$mean, 
                   test=uni_results$pvalue, 
                   theta=theta, method="BH")

plot(theta, rejBH, type="l",
     xlab=expression(theta), ylab="number of rejections")


theta = seq(from=0, to=0.99, by=0.02)
filterChoices = data.frame(
  `mean` = uni_results$mean,
  `geneID` = 1:nrow(uni_results),
  `sd` = uni_results$sd,
  `t.stat` = uni_results$test.stat
)
rejChoices = sapply(filterChoices, function(f)
  filtered_R(alpha=0.1, filter=f, test=uni_results$pvalue, theta=theta, method="BH"))

library("RColorBrewer")
myColours = brewer.pal(ncol(filterChoices), "Set1")
matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
        xlab=expression(theta), ylab="number of rejections")
legend("topleft", legend=colnames(filterChoices), fill=myColours)

rejChoices

theta = theta[which.max(rejChoices[,"mean"])]
pass = with(uni_results, mean > quantile(mean, theta))

h1 = hist(uni_results$pvalue[!pass], breaks=50, plot=FALSE)
h2 = hist(uni_results$pvalue[pass], breaks=50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


resFilt = uni_results[pass,]
orderInPlot = order(resFilt$pvalue)
showInPlot = (resFilt$pvalue[orderInPlot] <= 0.06)
alpha = 0.1

plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]),
     xlim=c(0,4000))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)

## ---- KS-filter ----


source("https://raw.githubusercontent.com/sahirbhatnagar/eclust/master/k_filter.R")

X <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ)][,brain_probes, with = F]
Y <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ)][,"WASI_Full_Scale_IQ", with = F]

colnames(DT_with_pheno)[1:30]

obj <- k.filter(x = X[,101:200, with=F],
                y = Y$WASI_Full_Scale_IQ, response.type = "continuous", method = "fused")

ks.test()
obj$k.rank
k.filter.single

warnings()

order(obj$k.rank)[1:5]

dat <- data.frame(x=t(DT.placenta.all[order(obj$k.rank)[1:5],,drop=F]), y=DT.pheno.placenta$`imc z-score`)

lm(y ~ ., dat) %>% summary()




## ---- data ----

X <- DT_with_pheno[!is.na(WASI_Full_Scale_IQ), c("Subject_ID","age_binary","Subject_Gender","Site_Location",brain_probes), with = F]

X %>% dim

filterd_probes <- uni_results[order(test.stat,decreasing = T)[1:3000]]$probe

X_all_reduced <- X[, c("Subject_ID","age_binary","Subject_Gender","Site_Location",
                       filterd_probes), with = F]

X_exposed_reduced <- X_all_reduced[age_binary %in% "(11.3,Inf]"]
X_unexposed_reduced <- X_all_reduced[age_binary %in% "(4.8,11.3]"]

TOM_all <- TOMsimilarityFromExpr(X_all_reduced[,filterd_probes,with=F], nThreads = 35)
TOM_exposed <- TOMsimilarityFromExpr(X_exposed_reduced[,filterd_probes,with=F], nThreads = 35)
TOM_unexposed <- TOMsimilarityFromExpr(X_unexposed_reduced[,filterd_probes,with=F], nThreads = 35)
TOM_diff <- abs(TOM_exposed-TOM_unexposed)

## ---- not-used ----
pheatmap::pheatmap(TOM_all, color = viridis(100), clustering_method = "average",
                   show_rownames = F, show_colnames = F)
