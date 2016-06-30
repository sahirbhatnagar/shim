## ---- required-packages ----

getPckg <- function(pckg) install.packages(pckg, repos = "http://cran.r-project.org")
getPckgBioc <- function(pckg) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(pckg)
}

pckg = try(require(magrittr))
if(!pckg) {
  cat("Installing 'magrittr' from CRAN\n")
  getPckg("magrittr")
  require(magrittr)
}

pckg = try(require(data.table))
if(!pckg) {
  cat("Installing 'data.table' from CRAN\n")
  getPckg("data.table")
  require(data.table)
}

pckg = try(require(plyr))
if(!pckg) {
  cat("Installing 'plyr' from CRAN\n")
  getPckg("plyr")
  require(plyr)
}

pckg = try(require(dplyr))
if(!pckg) {
  cat("Installing 'dplyr' from CRAN\n")
  getPckg("dplyr")
  require(dplyr)
}

pckg = try(require(WGCNA))
if(!pckg) {
  cat("Installing 'WGCNA' from CRAN\n")
  getPckg("WGCNA")
  require(WGCNA)
}
allowWGCNAThreads()

pckg = try(require(glmnet))
if(!pckg) {
  cat("Installing 'glmnet' from CRAN\n")
  getPckg("glmnet")
  require(glmnet)
}

pckg = try(require(gglasso))
if(!pckg) {
  cat("Installing 'gglasso' from CRAN\n")
  getPckg("gglasso")
  require(gglasso)
}

pckg = try(require(grpreg))
if(!pckg) {
  cat("Installing 'grpreg' from CRAN\n")
  getPckg("grpreg")
  require(grpreg)
}

pckg = try(require(dynamicTreeCut))
if(!pckg) {
  cat("Installing 'dynamicTreeCut' from CRAN\n")
  getPckg("dynamicTreeCut")
  require(dynamicTreeCut)
}

pckg = try(require(ncvreg))
if(!pckg) {
  cat("Installing 'ncvreg' from CRAN\n")
  getPckg("ncvreg")
  require(ncvreg)
}

pckg = try(require(PMA))
if(!pckg) {
  cat("Installing 'PMA' from CRAN\n")
  getPckg("PMA")
  require(PMA)
}

pckg = try(require(eclust))
if(!pckg) {
  cat("Installing 'eclust' from CRAN\n")
  getPckg("eclust")
  require(eclust)
}

pckg = try(require(doMC))
if(!pckg) {
  cat("Installing 'doMC' from CRAN\n")
  getPckg("doMC")
  require(doMC)
}

pckg = try(require(cluster))
if(!pckg) {
  cat("Installing 'cluster' from CRAN\n")
  getPckg("cluster")
  require(cluster)
}

pckg = try(require(Matrix))
if(!pckg) {
  cat("Installing 'Matrix' from CRAN\n")
  getPckg("Matrix")
  require(Matrix)
}

# pckg = try(require(protoclust))
# if(!pckg) {
#   cat("Installing 'protoclust' from CRAN\n")
#   getPckg("protoclust")
#   require(protoclust)
# }

pckg = try(require(stringr))
if(!pckg) {
  cat("Installing 'stringr' from CRAN\n")
  getPckg("stringr")
  require(stringr)
}

pckg = try(require(ggplot2))
if(!pckg) {
  cat("Installing 'ggplot2' from CRAN\n")
  getPckg("ggplot2")
  require(ggplot2)
}

pckg = try(require(latex2exp))
if(!pckg) {
  cat("Installing 'latex2exp' from CRAN\n")
  getPckg("latex2exp")
  require(latex2exp)
}

pckg = try(require(factoextra))
if(!pckg) {
  cat("Installing 'factoextra' from CRAN\n")
  getPckg("factoextra")
  require(factoextra)
}

pckg = try(require(FDb.InfiniumMethylation.hg19))
if(!pckg) {
  cat("Installing 'FDb.InfiniumMethylation.hg19' from Bioconductor\n")
  getPckgBioc("FDb.InfiniumMethylation.hg19")
  require(FDb.InfiniumMethylation.hg19)
}

pckg = try(require(genefilter))
if(!pckg) {
  cat("Installing 'genefilter' from Bioconductor\n")
  getPckgBioc("genefilter")
  require(genefilter)
}

# biocLite(c("qvalue","lumi","lumiHumanAll.db","lumiHumanIDMapping","illuminaHumanv4.db"))
# install.packages("bouchard_0.1.tar.gz", repos = NULL, type = "source")
require(bouchard)


pckg = try(require(readxl))
if(!pckg) {
  cat("Installing 'readxl' from CRAN\n")
  getPckg("readxl")
  require(readxl)
}

pckg = try(require(bit64))
if(!pckg) {
  cat("Installing 'bit64' from CRAN\n")
  getPckg("bit64")
  require(bit64)
}

pckg = try(require(pheatmap))
if(!pckg) {
  cat("Installing 'pheatmap' from CRAN\n")
  getPckg("pheatmap")
  require(pheatmap)
}

pckg = try(require(viridis))
if(!pckg) {
  cat("Installing 'viridis' from CRAN\n")
  getPckg("viridis")
  require(viridis)
}

pckg = try(require(R.matlab))
if(!pckg) {
  cat("Installing 'R.matlab' from CRAN\n")
  getPckg("R.matlab")
  require(R.matlab)
}

pckg = try(require(zoo))
if(!pckg) {
  cat("Installing 'zoo' from CRAN\n")
  getPckg("zoo")
  require(zoo)
}

pckg = try(require(RcppMLPACK))
if(!pckg) {
  cat("Installing 'RcppMLPACK' from CRAN\n")
  getPckg("RcppMLPACK")
  require(RcppMLPACK)
}

pckg = try(require(pls))
if(!pckg) {
  cat("Installing 'pls' from CRAN\n")
  getPckg("pls")
  require(pls)
}

