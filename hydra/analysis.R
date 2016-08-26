##################################
# R source code file for analysing eclust simulation study results on HYDRA cluster
# There is a separate analysis file for analyzing results from MAMMOUTH cluster
# Created by Sahir, May 11, 2016
# Updated: August 25, 2016
##################################
rm(list = ls())
source("~/git_repositories/eclust-simulation-aug2016/packages.R")
source("~/git_repositories/eclust-simulation-aug2016/functions.R")
options(digits = 2, scipen=999)

# note that some of the stability measures are NA in the Eclust, because both mclusters got 0 coefficients, therefore
# the intersect and union command give NAs
## ---- data ----

#source(paste(Sys.getenv("HOME"),"eclust/bin/simulation/sim_functions.R", sep = "/"))

col.names <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/colnames_nostab_nouni", header = F)$V1

DT <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/sim-hydra-results-p3000-no-stability", stringsAsFactors = FALSE) %>%
  setnames(col.names)

DT[, `:=`(simulation = 1:nrow(DT))]

# this still has all the raw data, but melted
options(warning.length = 8170)
DT.long <- DT %>%
  reshape2::melt(id.vars = c("simulation",colnames(DT)[1:20])) %>%
  tidyr::separate(variable, c("method", "summary", "model", "interaction", "measure"), convert = T) %>%
  as.data.table


DT.long$method %>% table
DT.long$model %>% table
DT.long$summary %>% table

levels <- c("pen", "clust", "Eclust")

DT.long[, `:=`(method = factor(method, levels = levels))]

DT.long[, table(measure)]
DT.long[measure=="FPR"][, hist(value)]
DT.long[measure=="TPR"][, hist(value)]


# this takes the mean by method across all simulations in a given method
DT.summary <- DT.long %>%
  tidyr::unite(name, summary, model) %>%
  summarySE(measurevar = "value", 
            groupvars = c("rho","p","SNR","n","nActive","Ecluster_distance","betaMean","alphaMean",
                          "measure","method","name"), 
            na.rm = TRUE) %>%
  as.data.table


DT.summary[, table(name)]
DT.summary[, table(rho)]
DT.summary[, table(p)]
DT.summary[, table(n)]
DT.summary[, table(nActive)]
DT.summary[, table(Ecluster_distance)]
DT.summary[, table(SNR)]
DT.summary[, table(betaMean)]
DT.summary[, table(alphaMean)]

levels.name <- c("na_elasticnet", "na_lasso","avg_elasticnet",
                 "avg_lasso", "avg_shim","pc_elasticnet",
                 "pc_lasso", "pc_shim")

labels.name <- c("elasticnet", "lasso","avg_elasticnet",
                 "avg_lasso", "avg_shim","pc_elasticnet",
                 "pc_lasso", "pc_shim")

DT.summary[,`:=`(name = factor(name, levels = levels.name, labels = labels.name))]

DT.summary$name %>% table
# TPR vs Shat group size vs rho -------------------------------------------------------------

## ---- tpr-vs-shat ----
p <- DT.long[measure %in% c("TPR","Shat")] %>%
  tidyr::spread(measure, value) %>%
  tidyr::unite(name, summary, model)
p[,`:=`(method=factor(method, levels = levels),name = factor(name, levels = levels.name, labels = labels.name) )]

p %>%
  ggplot(aes(x = Shat, y = TPR, color=method)) +
  geom_point(size=2.5, aes(shape=method)) +
  #geom_rug()+
  facet_grid(alphaMean+SNR~rho) +
  #facet_grid(name~size+rho)+
  #background_grid(major = "xy", minor = "xy")+
  theme_bw()+
  ylab("true positive rate")+
  xlab("number of non-zero estimated coefficients")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))+
  theme(legend.position = "top")
#ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/TPR_vs_Shat.png", sep = "/"))

# TPR vs FPR group size vs rho  -------------------------------------------------------------

## ---- tpr-vs-fpr ----
#SNR=0.2
p <- DT.long[measure %in% c("TPR","FPR")] %>%
  tidyr::spread(measure, value) %>%
  tidyr::unite(name, summary, model)
p[,`:=`(method=factor(method, levels = levels),name = factor(name, levels = levels.name, labels = labels.name) )]

p[name %in% c("lasso","elasticnet", "avg_elasticnet", "avg_shim")][SNR==0.2] %>%
  ggplot(aes(x = FPR, y = TPR, color=method)) +
  geom_jitter(size=2.5, aes(shape=name)) +
  #geom_rug()+
  facet_grid(alphaMean~rho) +
  #facet_grid(name~size+rho)+
  #background_grid(major = "xy", minor = "xy")+
  theme_bw()+
  ylab("true positive rate")+
  xlab("false positive rate")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))+
  theme(legend.position = "top")
# ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/TPR_vs_Shat.png", sep = "/"))
ggsave("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/tpr_fpr_snr02.png",
       width = 11, height = 8)

#SNR=1
p <- DT.long[measure %in% c("TPR","FPR")] %>%
  tidyr::spread(measure, value) %>%
  tidyr::unite(name, summary, model)
p[,`:=`(method=factor(method, levels = levels),name = factor(name, levels = levels.name, labels = labels.name) )]

p[name %in% c("lasso","elasticnet", "avg_elasticnet", "avg_shim")][SNR==1] %>%
  ggplot(aes(x = FPR, y = TPR, color=method)) +
  geom_jitter(size=2.5, aes(shape=name)) +
  #geom_rug()+
  facet_grid(alphaMean~rho) +
  #facet_grid(name~size+rho)+
  #background_grid(major = "xy", minor = "xy")+
  theme_bw()+
  ylab("true positive rate")+
  xlab("false positive rate")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))+
  theme(legend.position = "top")
# ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/TPR_vs_Shat.png", sep = "/"))
ggsave("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/tpr_fpr_snr1.png",
       width = 11, height = 8)

#SNR=2
p <- DT.long[measure %in% c("TPR","FPR")] %>%
  tidyr::spread(measure, value) %>%
  tidyr::unite(name, summary, model)
p[,`:=`(method=factor(method, levels = levels),name = factor(name, levels = levels.name, labels = labels.name) )]

p[name %in% c("lasso","elasticnet", "avg_elasticnet", "avg_shim")][SNR==2] %>%
  ggplot(aes(x = FPR, y = TPR, color=method)) +
  geom_jitter(size=2.5, aes(shape=name)) +
  #geom_rug()+
  facet_grid(alphaMean~rho) +
  #facet_grid(name~size+rho)+
  #background_grid(major = "xy", minor = "xy")+
  theme_bw()+
  ylab("true positive rate")+
  xlab("false positive rate")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))+
  theme(legend.position = "top")
# ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/TPR_vs_Shat.png", sep = "/"))
ggsave("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/tpr_fpr_snr2.png",
       width = 11, height = 8)


# MSE ---------------------------------------------------------------------
# pdf("protocol_simulation/plots/mse.pdf",width = 11, height = 8.5 )

## ---- mse ----
pd <- position_dodge(width = 0.7) # move them .05 to the left and right

ggplot(DT.summary[measure=="mse"][SNR==0.2],
  aes(x = name, y = mean, colour = method)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,size=1, position=pd) +
  geom_point(position=pd, size=3)+
  guides(color=guide_legend(title="method"))+
  xlab("model")+
  ylab("test set mean squared error")+
  #ylab(TeX("average MSE (1000 simulations)"))+
  facet_grid(alphaMean+SNR~rho, scales="free")+
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))
#ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/MSE_1_35.png", sep = "/"))
ggsave("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/mse_snr02.png",
       width = 11, height = 8)



ggplot(DT.summary[measure=="mse"][SNR==1],
       aes(x = name, y = mean, colour = method)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,size=1, position=pd) +
  geom_point(position=pd, size=3)+
  guides(color=guide_legend(title="method"))+
  xlab("model")+
  ylab("test set mean squared error")+
  #ylab(TeX("average MSE (1000 simulations)"))+
  facet_grid(alphaMean+SNR~rho, scales="free")+
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))
#ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/MSE_1_35.png", sep = "/"))
ggsave("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/mse_snr1.png",
       width = 11, height = 8)



ggplot(DT.summary[measure=="mse"][SNR==2],
       aes(x = name, y = mean, colour = method)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,size=1, position=pd) +
  geom_point(position=pd, size=3)+
  guides(color=guide_legend(title="method"))+
  xlab("model")+
  ylab("test set mean squared error")+
  #ylab(TeX("average MSE (1000 simulations)"))+
  facet_grid(alphaMean+SNR~rho, scales="free")+
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))
#ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/MSE_1_35.png", sep = "/"))
ggsave("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/mse_snr2.png",
       width = 11, height = 8)



# JACC, Spearman, Pearson (USED IN PROTOCOL)---------------------------------------------------------------------

# pdf("protocol_simulation/plots/jacc.pdf",width = 11, height = 8.5 )

## ---- jacc ----
#DT.summary[measure=="jacc"][is.na(mean), mean:=0]
#DT.summary[is.na(mean), mean:=0]
#DT.summary[is.na(mean)]

pd <- position_dodge(width = 0.7) # move them .05 to the left and right
ggplot(DT.summary[measure=="jacc"],
       aes(x = name, y = mean, colour = method)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,size=1, position=pd) +
  geom_point(position=pd, size=3)+
  guides(color=guide_legend(title="method"))+
  xlab("model")+
  ylab("Jaccard index")+
  #ylab(TeX("average MSE (1000 simulations)"))+
  facet_grid(p~rho, scales="fixed")+
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))
#ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/jaccard.png", sep = "/"))

# JACC, Spearman, Pearson (USED IN PROTOCOL)---------------------------------------------------------------------
# pdf("protocol_simulation/plots/spearman.pdf",width = 11, height = 8.5 )

## ---- spearman ----
pd <- position_dodge(width = 0.7) # move them .05 to the left and right
ggplot(DT.summary[measure=="spearman"],
       aes(x = name, y = mean, colour = method)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,size=1, position=pd) +
  geom_point(position=pd, size=3)+
  guides(color=guide_legend(title="method"))+
  xlab("model")+
  ylab("Spearman correlation")+
  #ylab(TeX("average MSE (1000 simulations)"))+
  facet_grid(p~rho, scales="fixed")+
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))
#ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/spearman.png", sep = "/"))

# JACC, Spearman, Pearson (USED IN PROTOCOL)---------------------------------------------------------------------
# pdf("protocol_simulation/plots/pearson.pdf",width = 11, height = 8.5 )

## ---- pearson ----
pd <- position_dodge(width = 0.7) # move them .05 to the left and right
ggplot(DT.summary[measure=="pearson"],
       aes(x = name, y = mean, colour = method)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,size=1, position=pd) +
  geom_point(position=pd, size=3)+
  guides(color=guide_legend(title="method"))+
  xlab("model")+
  ylab("Pearson correlation")+
  #ylab(TeX("average MSE (1000 simulations)"))+
  facet_grid(SNR~rho, scales="fixed")+
  theme_bw()+
  theme(legend.position = "top")+
  theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
        axis.text.y  = element_text(size=15),
        axis.title.x = element_text(face="bold", colour="#990000", size=15),
        axis.title.y = element_text(face="bold", colour="#990000", size=20),
        title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
        strip.text = element_text(size=20))+
  theme(legend.key.width=unit(1, "inches"))
#ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/pearson.png", sep = "/"))







