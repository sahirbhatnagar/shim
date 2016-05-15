#' ################################# 
#' R source code file used to get new
#' phenotype measures for analysis Created by Sahir, March 17, 2016 Updated: 
#' NOTE: We received more BMI measures from Luigi in March 2016. This script is
#' to update the phenotype information NOTE: For IDs 31,172,222 I manually
#' updated the DG case status, sex, age_gestationnel in the
#' CovariableModified.csv I took their information from 
#' "/home/data1/share/greenwood.group/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/bloodcord/COVARIABLE_DATA/Donnees_Pheno_Enfants_2014-12-09"
#' This is an updated version of previous files: 
#' /mnt/GREENWOOD_BACKUP/share/sy/bouchard/methylation/phenotype.R,
#' phenotype_placenta.R 
#' Updated: May 13, 2016.
#' This script is being called by data_cleaning.R
#' Now on git branch bouchard-analysis in eclust repo
#' #################################

#rm(list = ls())
# library(data.table)
# library(readxl)

# get phenotypes for placenta  ---------------------------------------------------------

# get gestational diabetes phenotype..NOTE: the phenotype are the same as bloodcoord, but the sample ID names
# are different which is why you need to load this one for the code to work
#DT.pheno.placenta <- fread("/mnt/GREENWOOD_BACKUP/share/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/placenta/COVARIABLE_DATA/CovariableModified.csv")
DT.pheno.placenta <- fread("~/Documents/bouchard_data/CovariableModified.csv")
set(DT.pheno.placenta,i = NULL, j = "ID", value = paste(DT.pheno.placenta[["Sentrix_ID"]], DT.pheno.placenta[["Sentrix_Position"]], sep = "_"))

setkey(DT.pheno.placenta,ID)
#remove extra 61F_Rep
DT.pheno.placenta <- DT.pheno.placenta[-which(ID=="5975819046_R06C02")]

nrow(DT.pheno.placenta)

#setnames(DT.pheno.placenta, c("Diabetes", "Age.gestationnel"), c("diabetes","age"))
DT.pheno.placenta[,case := factor(Diabete_Gest, levels =c ("NGT","DG"))]
DT.pheno.placenta[,Sexe := factor(Sexe, levels = c("f","m"))]

DT.pheno.placenta <- DT.pheno.placenta[,c("Sample_ID", "ID","Diabete_Gest","case","Age_gestationnel","AgeMois","BMI","ZScoreBMI",
                                          "percentFAT","Tricep","Bicep","Sous_Scapulaire","Iliaque","Sexe"), with=F]

DT.pheno.placenta <- DT.pheno.placenta[!is.na(case)]

nrow(DT.pheno.placenta)

# this checks that all non-missing percentFAT ID's are found in the methylation data
#DT.pheno.placenta[!is.na(percentFAT)]$ID %in% colnames(DT)

DT.pheno.placenta[, Sample_ID_num := as.double(gsub("F","",Sample_ID))]
setkey(DT.pheno.placenta, Sample_ID_num)
key(DT.pheno.placenta)

# get phenotypes for cordblood (not used as of May 2016, I am focusing on placenta) ---------------------------------------------------------

# DT.pheno.cordblood <- fread("/mnt/GREENWOOD_BACKUP/share/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/bloodcord/COVARIABLE_DATA/CovariableModified.csv")
# set(DT.pheno.cordblood, i = NULL, j = "ID", 
#     value = paste(DT.pheno.cordblood[["Sentrix_ID"]],DT.pheno.cordblood[["Sentrix_Position"]],sep="_"))
# 
# # as per Andre Anne's email, ID 59,95,155 do not have cordblood methylation data
# # note also that sample 222 did not pass QC filter results as per
# # setwd("/home/data1/share/greenwood.group/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/bloodcord/QC_REPORT/excluded_samples.R")
# # so ID 6229050136_R02C01 shows up in DT.pheno.cordblood2 but not in the raw methylation data
# 
# setkey(DT.pheno.cordblood,ID)
# #remove extra 61F_Rep
# DT.pheno.cordblood <- DT.pheno.cordblood[-which(ID=="6229050136_R06C02")]
# 
# #setnames(DT.pheno.cordblood, c("Diabetes", "Age.gestationnel"), c("diabetes","age"))
# DT.pheno.cordblood[,case := factor(Diabete_Gest, levels=c("NGT","DG"))]
# DT.pheno.cordblood[,Sexe := factor(Sexe, levels=c("f","m"))]
# 
# DT.pheno.cordblood <- DT.pheno.cordblood[,c("Sample_ID", "ID","Diabete_Gest","case","Age_gestationnel","AgeMois","BMI","ZScoreBMI",
#                                             "percentFAT","Tricep","Bicep","Sous_Scapulaire","Iliaque","Sexe"), with=F]
# 
# # as per Andre Anne's email, ID 59,95,155 do not have cordblood methylation data
# # note also that sample 222 did not pass QC filter results as per
# # setwd("/home/data1/share/greenwood.group/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/bloodcord/QC_REPORT/excluded_samples.R")
# # so ID 6229050136_R02C01 shows up in DT.pheno2 but not in the raw methylation data
# # remove blank Sample_ID and ID 222F 
# DT.pheno.cordblood <- DT.pheno.cordblood[!is.na(case)][Sample_ID!="222F"]
# 
# DT.pheno.cordblood[, Sample_ID_num := as.double(gsub("F","",Sample_ID))]
# setkey(DT.pheno.cordblood, Sample_ID_num)
# key(DT.pheno.cordblood)

# get newer phenotype data sent in March 2016 -----------------------------

DT.pheno.placenta[is.na(ZScoreBMI)]
# DT.pheno.cordblood[is.na(ZScoreBMI)]

# DT.pheno.new <- as.data.table(readxl::read_excel("/mnt/GREENWOOD_BACKUP/share/PROJECTS/Luigi_Bouchard/Methylation/10-06-2014/placenta/COVARIABLE_DATA/phénotypes enfants 5 ans ZPREG.xlsx"))
DT.pheno.new <- as.data.table(readxl::read_excel("~/Documents/bouchard_data/phenotypes enfants 5 ans ZPREG.xlsx"))
DT.pheno.new <- DT.pheno.new[!is.na(ID_zpreg)]
DT.pheno.new[,str(ID_zpreg)]
setkey(DT.pheno.new, ID_zpreg)

# this shows that the points match up in each dataset
DT.pheno.new[DT.pheno.placenta][, plot(`imc z-score`, ZScoreBMI)]
# DT.pheno.new[DT.pheno.cordblood][, points(`imc z-score`, ZScoreBMI, pch=19)]
abline(a = 0, b = 1)


DT.pheno.new[DT.pheno.placenta][is.na(`imc z-score`)]
DT.pheno.new[DT.pheno.placenta][is.na(ZScoreBMI)]

# DT.pheno.new[DT.pheno.cordblood][is.na(`imc z-score`)]
# DT.pheno.new[DT.pheno.cordblood][is.na(ZScoreBMI)]


DT.pheno.placenta <- DT.pheno.new[DT.pheno.placenta]
# DT.pheno.cordblood <- DT.pheno.new[DT.pheno.cordblood]
# key(DT.pheno.cordblood)
key(DT.pheno.placenta)

# Load covariate data -----------------------
# we didnt have the right phenotype file originally, so we had to map three files in order to get the 
# proper phenotypes with the ID's that occur in the expression data
#setwd("~/share/greenwood.group/PROJECTS/Luigi_Bouchard/Expression")

# DT.cov <- fread("/mnt/GREENWOOD_BACKUP/share/PROJECTS/Luigi_Bouchard/Expression/WholeGenomeSampleSheet lBouchard001modified.csv")
DT.cov <- fread("~/Documents/bouchard_data/WholeGenomeSampleSheet lBouchard001modified.csv")
DT.cov[,IDNEW := paste(Sentrix_ID,Sentrix_Position,sep="_")]
setkey(DT.cov, Sample_Name)

# DT.id <- fread("/mnt/GREENWOOD_BACKUP/share/PROJECTS/Luigi_Bouchard/Expression/IDZPREG_EXPRESSION.csv")
DT.id <- fread("~/Documents/bouchard_data/IDZPREG_EXPRESSION.csv")
setkey(DT.id,IDZPREG)

DT.temp <- DT.id[DT.pheno.placenta]
setkey(DT.temp, ID)

# DT.temp2 <- DT.id[DT.pheno.cordblood]
# setkey(DT.temp2, ID)

# final phenotype dataset ------------------------------------------------
DT.pheno.placenta <- DT.cov[DT.temp]
# DT.pheno.cordblood <- DT.cov[DT.temp2]


#write.table(DT.pheno.expression,"~/share/sy/bouchard/expression/probes/DT.pheno.expression",row.names=FALSE)
#DT.pheno.expression
rm(DT.cov,DT.id,DT.temp, DT.pheno.new)

colnames(DT.pheno.placenta)
DT.pheno.placenta <- DT.pheno.placenta[,c("Sample_Name","IDNEW", "IDZPREG","age (mois)", "poids z-score", "imc", "imc z-score", "commentaire",
                                          "Sample_ID","i.ID","case", "Age_gestationnel", "Sexe"), with=F]

# colnames(DT.pheno.cordblood)
# DT.pheno.cordblood <- DT.pheno.cordblood[,c("Sample_Name","IDNEW", "IDZPREG","age (mois)", "poids z-score", "imc", "imc z-score", "commentaire",
#                                             "Sample_ID","i.ID","case", "Age_gestationnel", "Sexe"), with=F]








