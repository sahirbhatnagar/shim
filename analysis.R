##################################
# R source code file for data manipulation. This script
# imports the expr. and meth data along with annotations
# run filter_cca.R after running this code to filter the expression
# and methylation data and get the canonical variables (X_1 w_1 ...)
# canonical vectors are the weights (w_1...)
# Created by Sahir, March 21
# Updated:
# NOTE: This directory is more upto data than bouchard/scripts/ccareport
# I am now using version control in bouchard/cca/
# because it was getting too messy in ccareport
##################################

rm(list=ls())
source("functions.R")
source("packages.R")
source("data_cleaning.R")


## ---- corScor ----
require(doMC)
registerDoMC(cores = 3)
cor_scor_placenta <- bigcorPar(data.all = t(DT.placenta.all[1:10000,]),
                               data.e0 = t(DT.placenta.ngd[1:10000,]),
                               data.e1 = t(DT.placenta.gd[1:10000,]),
                               alpha = 2, threshold = 1, nblocks = 100,
                               ncore = 3)

cor_scor_placenta$score %>% as.numeric %>% hist

probes_placenta <- unique(c(cor_scor_placenta$gene1, cor_scor_placenta$gene2))


## ---- list-of-cg-sites ----


write.table(DT.placenta[probes_placenta,c("rn","nearestGeneSymbol","CHR"), with=F],file="probes_placenta.txt", quote = F, row.names = F)
colnames(DT.placenta)


cor_scor <- as.data.table(cor_scor_placenta)
cor_scor_m <- melt(cor_scor, measure.vars = c("gene1","gene2"))
setkey(cor_scor_m, "value")
unique(cor_scor_m)
cor_scor_unique <- unique(cor_scor_m)
cor_scor_unique %>% str
cor_scor_unique[,score := as.numeric(score)]


t1 <- DT.placenta[probes_placenta,c("rn","nearestGeneSymbol","CHR"), with=F]
setkey(t1, "rn")

t_final <- t1[cor_scor_unique[, c("score","value"), with = F]]

write.table(t_final,file="probes_placenta_with_score.txt", quote = F, row.names = F)

as.data.table(DT.placenta.all)

write.table(DT.placenta[,c("rn","nearestGeneSymbol","CHR" ), with=F],
            file = "probes_placenta_used_in_correlations.txt",
            quote=F, row.names = F)

## ---- rand-Index ----


p = length(probes_placenta)
res <- vector("list",3)
cor.matrix.x0 <- corFast(t(DT.placenta.ngd[probes_placenta,]))
cor.matrix.x1 <- corFast(t(DT.placenta.gd[probes_placenta,]))
corr <- corFast(t(DT.placenta.all[probes_placenta,]))

res <- lapply(c("corr", "corr0","corr1"), function(i) {
  
  distance <- switch(i,
                     corr = corr,
                     corr0 = cor.matrix.x0,
                     corr1 = cor.matrix.x1,
                     diffcorr = diff.corr,
                     difftom = diff.tom,
                     tom0 = tom.matrix.x0,
                     tom1 = tom.matrix.x1,
                     tom = tom)
  
  cl <- hclust(as.dist(if (i %in% c("diffcorr","difftom")) distance else 1-distance), method = "ward.D2")
  cuttree <- dynamicTreeCut::cutreeDynamic(cl,
                                           distM = if (i %in% c("diffcorr","difftom")) distance else 1-distance,
                                           deepSplit = 1)
  # check if all cluster groups are 0 which means no cluster assignment and everyone is in their own group
  clusters <- data.table(cpg = colnames(distance), cluster = if (all(cuttree==0)) 1:p else cuttree)
  setkey(clusters,"cpg")
  setnames(clusters, "cluster", i)
  clusters
})

Venn <- res[[1]][res[[2]]][res[[3]]]


library(mclust)

adjustedRandIndex(Venn$corr, Venn$corr0)
adjustedRandIndex(Venn$corr, Venn$corr1)
adjustedRandIndex(Venn$corr1, Venn$corr0)





## ---- eclust ----

sigProbes <- read.table("/mnt/GREENWOOD_BACKUP/share/sy/interactions/coexpression/bouchard/probes_placenta_with_score.txt",
                        header = T, stringsAsFactors = FALSE)

DT <- DT.placenta.all[sigProbes$rn[sigProbes$rn %in% rownames(DT.placenta.all)],] %>%
  
  sigProbes$rn[sigProbes$rn %in% rownames(DT.placenta.all)]