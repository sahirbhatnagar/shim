## ---- required-packages ----

getPckg <- function(pckg) install.packages(pckg, repos = "http://cran.r-project.org")
getPckgBioc <- function(pckg) biocLite(pkgs = pckg)


pckg = try(require(devtools))
if(!pckg) {
  cat("Installing 'devtools' from CRAN\n")
  getPckg("devtools")
  require(devtools)
}

# install_github('sahirbhatnagar/eclust')


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

source("http://bioconductor.org/biocLite.R")

pckg = try(require(Biobase))
if(!pckg) {
  cat("Installing 'Biobase' from Bioconductor\n")
  getPckgBioc("Biobase")
  require(Biobase)
}

pckg = try(require(IRanges))
if(!pckg) {
  cat("Installing 'IRanges' from Bioconductor\n")
  getPckgBioc("IRanges")
  require(IRanges)
}

pckg = try(require(AnnotationDbi))
if(!pckg) {
  cat("Installing 'AnnotationDbi' from Bioconductor\n")
  getPckgBioc("AnnotationDbi")
  require(AnnotationDbi)
}

pckg = try(require(impute))
if(!pckg) {
  cat("Installing 'impute' from Bioconductor\n")
  getPckgBioc("impute")
  require(impute)
}

pckg = try(require(preprocessCore))
if(!pckg) {
  cat("Installing 'preprocessCore' from Bioconductor\n")
  getPckgBioc("preprocessCore")
  require(preprocessCore)
}

pckg = try(require(GO.db))
if(!pckg) {
  cat("Installing 'GO.db' from Bioconductor\n")
  getPckgBioc("GO.db")
  require(GO.db)
}

pckg = try(require(WGCNA))
if(!pckg) {
  cat("Installing 'WGCNA' from CRAN\n")
  getPckg("WGCNA")
  require(WGCNA)
}

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

pckg = try(require(viridis))
if(!pckg) {
  cat("Installing 'viridis' from CRAN\n")
  getPckg("viridis")
  require(viridis)
}

pckg = try(require(pheatmap))
if(!pckg) {
  cat("Installing 'pheatmap' from CRAN\n")
  getPckg("pheatmap")
  require(pheatmap)
}

pckg = try(require(DT))
if(!pckg) {
  cat("Installing 'DT' from CRAN\n")
  getPckg("DT")
  require(DT)
}

pckg = try(require(caret))
if(!pckg) {
  cat("Installing 'caret' from CRAN\n")
  getPckg("caret")
  require(caret)
}

pckg = try(require(earth))
if(!pckg) {
  cat("Installing 'earth' from CRAN\n")
  getPckg("earth")
  require(earth)
}


pckg = try(require(cowplot))
if(!pckg) {
  cat("Installing 'cowplot' from CRAN\n")
  getPckg("cowplot")
  require(cowplot)
}

# pckg = try(require(fields))
# if(!pckg) {
#   cat("Installing 'fields' from CRAN\n")
#   getPckg("fields")
#   require(fields)
# }

pckg = try(require(RColorBrewer))
if(!pckg) {
  cat("Installing 'RColorBrewer' from CRAN\n")
  getPckg("RColorBrewer")
  require(RColorBrewer)
}

pckg = try(require(pROC))
if(!pckg) {
  cat("Installing 'pROC' from CRAN\n")
  getPckg("pROC")
  require(pROC)
}
