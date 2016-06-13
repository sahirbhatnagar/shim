##################################
# R source code file for importing and cleaning nihpd dataset
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


# DT <- R.matlab::readMat("~/Dropbox/PhD/data/NIHPD/NIHP_Cortical_Thickness.mat")
dat <- R.matlab::readMat("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_data/NIHP_Cortical_Thickness.mat")

str(dat)

DT <- dat$Matrix.All.81924 %>% 
  as.data.table() %>% 
  setnames(c("V1","V2","V3","V4"), c("ID", "age","F","M")) %>% 
  setkey(ID)
DT[,`:=`(age=age/365.25)]
# categorize the age as per Cereb.Cortex paper by Budha
DT[, fivenum(age)]
DT[, age_binary := cut(age, c(4.8, 11.3, Inf))]
DT[, table(age_binary)]
DT[,`:=`(freq=.N), by = ID]
DT[, table(freq, useNA = "always")]
DT[,1:4, with=F]

brain_probes <- grep("V\\d*", colnames(DT),  value = TRUE)
# brain_probes
# length(brain_probes)
# setdiff(colnames(DT), brain_probes)

setcolorder(DT, c("ID","age","age_binary","F","M","freq", brain_probes))

# data for the first timepoint (this means that its the first obseration for that person)
DT1 <- DT[DT[unique(DT), , mult = "first", which = TRUE]]

# nrow(DT1)
# dim(DT1)
# DT1[, table(ID)]
# DT1[, unique(ID)] %>% length
# DT1[,1:7, with=F]

# rm(pheno,pheno1)
# pheno <- as.data.table(readxl::read_excel("~/Dropbox/PhD/data/NIHPD/NIHPD_Demographics_IQ.xls"))
# there are some subjects in this NIHPD_Demographics_IQ.xls file that are not in the NIHP_Cortical_Thickness.mat file
# dont use the age from here, use the age found in the  NIHP_Cortical_Thickness.mat file
pheno <- as.data.table(readxl::read_excel("/mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/eclust/rda/nihpd_data/NIHPD_Demographics_IQ.xls"))
setkey(pheno, Subject_ID)
data.table::setnames(pheno, "QC_050514: Scores (0 = Failed; 1 = Check; 2 = Passed)", "QC")
data.table::setnames(pheno, "Age_Days_Date_of_Visit_to_DOB", "age")
pheno[, age:=NULL]
pheno <- pheno[!is.na(Timepoint_Label)]
pheno[, time:=NULL]
pheno[, time:=vector("integer", nrow(pheno))]
pheno[, time := if(Timepoint_Label %in% c("v1","V1")) 1L 
      else if (Timepoint_Label %in% c("v2","V2")) 2L 
      else if (Timepoint_Label %in% c("v3","V3","v4")) 3L, by=list(Subject_ID, Timepoint_Label)]

# pheno[Subject_ID %in% pheno[Timepoint_Label %in% "v4"]$Subject_ID]
# pheno[Timepoint_Label %in% "v4"]$Subject_ID %in% DT1[,1:4,with=F]$ID

# pheno[, xtabs(~Subject_ID+time)]
# pheno[, xtabs(~time+Timepoint_Label)]

pheno[, Scanner_Serial_T1:=as.numeric(Scanner_Serial_T1)]
pheno[, `:=`(WASI_Full_Scale_IQ=as.numeric(WASI_Full_Scale_IQ), 
             WASI_Performance_IQ=as.numeric(WASI_Performance_IQ),
             WASI_Verbal_IQ=as.numeric(WASI_Verbal_IQ),
             Household_Income_Level=as.numeric(Household_Income_Level))]

# income is in the dataset at time 1, but sometimes at time 2 or 3, 
# so need to use the na.locf to carry obseervations forwards or backward
pheno[, income_binary:=NULL]
pheno[, table(HUD_Adjusted_Family_Income, useNA = "always")]
pheno[HUD_Adjusted_Family_Income=="."]
pheno[HUD_Adjusted_Family_Income==".", HUD_Adjusted_Family_Income := NA]

pheno[, table(HUD_Adjusted_Family_Income, useNA = "always")]
pheno[, HUD_Adjusted_Family_Income := na.locf(HUD_Adjusted_Family_Income, na.rm = FALSE, fromLast = FALSE), by = Subject_ID]
pheno[, HUD_Adjusted_Family_Income := na.locf(HUD_Adjusted_Family_Income, na.rm = FALSE, fromLast = TRUE), by = Subject_ID]
pheno[, table(HUD_Adjusted_Family_Income, useNA = "always")]
pheno[, table(Subject_ID, HUD_Adjusted_Family_Income, useNA="always")]

pheno[, income := factor(HUD_Adjusted_Family_Income, levels = c("Low", "Medium", "High")) ]
pheno[income %in% c("Low","Medium"), income_binary := "Low" ]
pheno[income %in% "High", income_binary := "High"]
pheno[, income_binary := factor(income_binary, levels = c("Low", "High"))]
pheno[, table(income_binary, useNA = "always")]

# also try to create binary income variable based on Household_Income_Level variable
# low: score from 1-7, high: score from 8-10
# this leads to perfectly balanced groups.. 164 in the low, 165 in the high

pheno[, table(Household_Income_Level, useNA="always")]
pheno[, Household_Income_Level := na.locf(Household_Income_Level, na.rm = FALSE, fromLast = FALSE), by = Subject_ID]
pheno[, Household_Income_Level := na.locf(Household_Income_Level, na.rm = FALSE, fromLast = TRUE), by = Subject_ID]

with(pheno,xtabs(~Household_Income_Level+income_binary))
pheno[Household_Income_Level %in% 1:7, income_binary2 := "Low"]
pheno[Household_Income_Level %in% 8:10, income_binary2 := "High"]
pheno[, income_binary2 := factor(income_binary2, levels = c("Low","High"))]
pheno[, xtabs(~Household_Income_Level + income_binary2)]

# pheno[, xtabs(~QC+time)]
# pheno[, table(time)]

# data for the first timepoint of each subject
pheno1 <- pheno[pheno[unique(pheno), , mult = "first", which = TRUE]]
pheno1[, table(income_binary, useNA = "always")]

# ll <- DT[,1:4,with=F]

# 57 subjects in the phenotype file who are not in the brain data file
# unique(pheno$Subject_ID)[which(unique(pheno$Subject_ID) %ni% ll$ID)]
# 
# all(dat$subj[,1] %in% ll$ID)
# any(dat$subj[,1] %in% 1009)

# pheno1[!is.na(WASI_Full_Scale_IQ) & QC > 0]
# pheno1[is.na(income_binary)]
# pheno1[!is.na(income_binary) & QC > 0]
# pheno1[!is.na(income_binary2) & QC >= 0]
# pheno1[!is.na(WASI_Full_Scale_IQ) & QC >= 0]
# pheno1[!is.na(WASI_Performance_IQ) & QC >= 0]
# pheno1[!is.na(WASI_Verbal_IQ) & QC >= 0]
# pheno1[!is.na(WASI_Verbal_IQ) & QC >= 0, hist(WASI_Verbal_IQ)]
# pheno1[!is.na(WASI_Verbal_IQ) & QC >= 0, hist(WASI_Full_Scale_IQ)]
# pheno1[!is.na(WASI_Verbal_IQ) & QC >= 0, hist(WASI_Performance_IQ)]
# pheno1[, .N, by = list(Subject_ID)][, table(N)]

# DT1[,1:4, with=F][pheno1[QC>=0]][, table(income_binary, useNA = "always")]
# DT1[,1:4, with=F][pheno1[QC>0]][, table(age_binary, useNA = "always")]
# pheno1[DT1[,1:4, with=F]][, table(age_binary, useNA = "always")]
# pheno1[DT1[,1:4, with=F]][, table(income_binary, useNA = "always")]
# 
# pheno1[DT1[,1:4, with=F]][!is.na(WASI_Full_Scale_IQ)]
# pheno1[DT1[,1:4, with=F]][!is.na(WASI_Performance_IQ)]
# pheno1[DT1[,1:4, with=F]][!is.na(WASI_Performance_IQ)]
# pheno1[DT1[,1:4, with=F]][!is.na(age_binary)]
# pheno[DT1[,1:4, with=F]][!is.na(income_binary2)]

# rm(DT_with_pheno)
DT_with_pheno <- pheno1[DT1]
dim(DT_with_pheno)

setcolorder(DT_with_pheno, c(setdiff(colnames(DT_with_pheno), brain_probes), brain_probes))

DT_with_pheno[, Subject_Gender := factor(Subject_Gender, levels = c("Male","Female"))]
DT_with_pheno[, Site_Location := factor(Site_Location)]
DT_with_pheno[, "Site_Location", with = F] %>% str

# convert age_binary into a numeric so that the fitting functions work
DT_with_pheno[, E := as.numeric(age_binary)-1]
DT_with_pheno[, table(E, age_binary)]
# DT_with_pheno[, age_binary ]
# DT_with_pheno[, table(income_binary, useNA = "always")]
# DT_with_pheno[, table(income_binary2, useNA = "always")]
# DT_with_pheno[!is.na(WASI_Full_Scale_IQ), .N]
# DT_with_pheno[!is.na(age_binary)] %>% dim
# 
# pheno1[, table(income_binary, useNA = "always")]

# this shows that there are some subject in pheno that are not in the data file
# pheno1$Subject_ID %in% DT[,1:5,with=F]$ID



# xyplot(QC ~ time , type = "l",  groups = Subject_ID, data=pheno)
# xyplot(Scanner_Serial_T1 ~ time , type = "l",  groups = Subject_ID, data=pheno)
# xyplot(age ~ time , type = "l",  groups = Subject_ID, data=pheno)
# xyplot(income ~ time , type = "l",  groups = Subject_ID, data=pheno)
# xyplot(WASI_Full_Scale_IQ ~ time , type = "l",  groups = Subject_ID, 
#        data=pheno[Subject_ID %in% sample(Subject_ID,10)])
# xyplot(WASI_Full_Scale_IQ ~ time , type = "l",  groups = Subject_ID, 
#        data=pheno)


