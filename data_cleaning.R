##################################
# R source code file for data manipulation of Bouchard Data
# There's expression data for placenta data only
# Created by Sahir, May 13, 2016
# Updated:
##################################


## ---- data-placenta ----

# source("/mnt/GREENWOOD_BACKUP/share/sy/bouchard/cca/scripts/multiplot.R")
# source("/mnt/GREENWOOD_BACKUP/share/sy/interactions/coexpression/bouchard/bigcorPar.R")

hm450 = FDb.InfiniumMethylation.hg19::get450k()

# load("/mnt/GREENWOOD_BACKUP/share/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/bloodcord/FILTERED_DATA/matrix_filter1.RData")
# load("~/Documents/methylation_bouchard/bloodcord/matrix_filter1.RData")
# DT.raw.cord <- as.data.table(filtered_matrix, keep.rownames = TRUE)
# setkey(DT.raw.cord,rn)
# n <- nrow(DT.raw.cord)
# # remove the following sample, as it is blank in covariate file
# DT.raw.cord[,"6229050123_R01C02" := NULL,]
# # remove the following sample as it is 61_F_Replicate
# DT.raw.cord[, "6229050136_R06C02" := NULL,]

# load("/mnt/GREENWOOD_BACKUP/share/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/placenta/FILTERED_DATA/matrix_filter1.RData")
load("~/Documents/bouchard_data/methylation/placenta/matrix_filter1.RData")

DT.raw.placenta <- as.data.table(filtered_matrix, keep.rownames = TRUE)
# Remove 61F_Rep replicated subject
DT.raw.placenta[,"5975819046_R06C02" := NULL, ]
setkey(DT.raw.placenta,rn)
n <- nrow(DT.raw.placenta)
DT.raw.placenta <- DT.raw.placenta[, lapply(.SD, function(x){(x*(n - 1) + 0.5) / n}), by = rn]

# DT.raw.cord[, "mean_methylation" := rowMeans(DT.raw.cord[,-1, with = FALSE])]
DT.raw.placenta[, "mean_methylation" := rowMeans(DT.raw.placenta[,-1, with = FALSE])]
# DT.raw.cord[, "sd_methylation" := rowSds(DT.raw.cord[,which(colnames(DT.raw.cord) %ni% c("rn", "mean_methylation")), with = FALSE])]
DT.raw.placenta[, "sd_methylation" := rowSds(DT.raw.placenta[,which(colnames(DT.raw.placenta) %ni% c("rn", "mean_methylation")), with = FALSE])]

# annotated dataset
DT.placenta <- cg.annotate(DT.raw.placenta[mean_methylation>=0.1 & mean_methylation <= 0.9])

# DT.cord <- cg.annotate(DT.raw.cord[mean_methylation>=0.1 & mean_methylation <= 0.9])

# # for celia
# write.table(DT.raw.placenta$sd_methylation, file = "sd_placenta.txt", quote = F, row.names = F, col.names = F)
# write.table(DT.raw.cord$sd_methylation, file = "sd_cord.txt", quote = F, row.names = F, col.names = F)
# g <- fread("sd_cord.txt")

# load phenotype
source("phenotype.R")

# ID of GD cases
GD <- DT.pheno.placenta[!is.na(`imc z-score`)][case == "DG"]$i.ID

# ID of non-GD cases
NGD <- DT.pheno.placenta[!is.na(`imc z-score`)][case == "NGT"]$i.ID


GD %in% colnames(DT.placenta) %>% sum
NGD %in% colnames(DT.placenta) %>% sum


str(DT.pheno.placenta)


#DT.placenta.gd <- DT.placenta[sd_methylation > quantile(DT.raw.placenta$sd_methylation,probs = c(0.4))][,colnames(DT.placenta) %in% GD, with = F] %>% as.matrix()

DT.placenta.gd <- DT.placenta[,colnames(DT.placenta) %in% GD, with = F] %>% as.matrix()
dim(DT.placenta.gd)
dimnames(DT.placenta.gd)[[1]] <- DT.placenta$rn

DT.placenta.ngd <- DT.placenta[,colnames(DT.placenta) %in% NGD, with = F] %>% as.matrix()
dim(DT.placenta.ngd)
dimnames(DT.placenta.ngd)[[1]] <- DT.placenta$rn

DT.placenta.all <- DT.placenta[,colnames(DT.placenta) %in% c(GD,NGD), with = F] %>% as.matrix()
dim(DT.placenta.all)
dimnames(DT.placenta.all)[[1]] <- DT.placenta$rn

rm(filtered_matrix, DT.raw.placenta)

setkey(DT.pheno.placenta, NULL)
setkey(DT.pheno.placenta, i.ID)

dim(DT.placenta)

# only keep the phenotypes for who we have methylation data
DT.pheno.placenta <- DT.pheno.placenta[match(dimnames(DT.placenta.all)[[2]],DT.pheno.placenta$i.ID)]


DT.placenta.all_datatable <- as.data.table(DT.placenta.all)
# reorder columns of methylation data to match the ordering of the phenotype data
data.table::setcolorder(DT.placenta.all_datatable, DT.pheno.placenta$i.ID)
dim(DT.placenta.all_datatable)

# check that ordering is correct
colnames(DT.placenta.all_datatable)==DT.pheno.placenta$i.ID

# convert back to matrix
DT.placenta.all <- as.matrix(DT.placenta.all_datatable)

