##################################
# R source code file for importing nihpd dataset
# 
# 
# 
# 
# Created by Sahir, May 16
# Updated:
# NOTE: 
# 
# 
##################################

source("packages.R")


DT <- R.matlab::readMat("~/Dropbox/PhD/data/NIHPD/NIHP_Cortical_Thickness.mat")


pheno <- as.data.table(readxl::read_excel("~/Dropbox/PhD/data/NIHPD/NIHPD_Demographics_IQ.xls"))
setkey(pheno, Subject_ID)



subID <- as.data.table(DT$subj)
setkey(subID, V1)

# the index of the first occurence for each subject
idIndex <- subID[unique(subID), , mult = "first", which = TRUE]

# categorize the age as per Cereb.Cortex paper by Budha
ageGroup <- cut(DT$age[idIndex,]/365.25,c(4.8, 11.3, Inf))


# This contains SubID, age in days, Female, Male, 81924 measurements
DT$Matrix.All.81924[,1:4]

DT$gen
DT$