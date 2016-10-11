##################################
# R source code file for analysing eclust simulation study results on HYDRA cluster
# There is a separate analysis file for analyzing results from MAMMOUTH cluster
# Created by Sahir, May 11, 2016
# Updated: September 11, 2016
##################################

rm(list = ls())
source("~/git_repositories/eclust-simulation-aug2016/packages.R")
source("~/git_repositories/eclust-simulation-aug2016/functions.R")
options(digits = 2, scipen=999)

# note that some of the stability measures are NA in the Eclust, because both mclusters got 0 coefficients, therefore
# the intersect and union command give NAs
## ---- data ----

col.names <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/colnames_stab_hydra-sim3-sept27.txt", header = F)$V1

DT <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/sim3-hydra-results-p5000-sept27", stringsAsFactors = FALSE) %>%
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
DT.long$cluster_distance %>% table
DT.long$Ecluster_distance %>% table


DT.long[, `:=`(method = factor(method, levels = c("mars", "clust", "Eclust"), labels = c("SEPARATE", "CLUST", "ECLUST")))]
DT.long[, `:=`(cluster_distance = factor(cluster_distance, levels = c("corr","tom"), labels = c("Correlation", "TOM")))]

DT.long[, table(measure)]
DT.long$cluster_distance %>% table
DT.long[, table(method)]
DT.long[measure=="FPR"][, hist(value)]
DT.long[measure=="TPR"][, hist(value)]


# this takes the mean by method across all simulations in a given method
# DT.summary <- DT.long %>%
#   tidyr::unite(name, summary, model) %>%
#   summarySE(measurevar = "value", 
#             groupvars = c("rho","p","SNR","n","nActive","Ecluster_distance","betaMean","alphaMean",
#                           "measure","method","name"), 
#             na.rm = TRUE) %>%
#   as.data.table
# 
# 
# DT.summary[, table(name)]
# DT.summary[, table(rho)]
# DT.summary[, table(p)]
# DT.summary[, table(n)]
# DT.summary[, table(nActive)]
# DT.summary[, table(Ecluster_distance)]
# DT.summary[, table(SNR)]
# DT.summary[, table(betaMean)]
# DT.summary[, table(alphaMean)]
# 
# levels.name <- c("na_elasticnet", "na_lasso","avg_elasticnet",
#                  "avg_lasso", "avg_shim","pc_elasticnet",
#                  "pc_lasso", "pc_shim")
# 
# labels.name <- c("elasticnet", "lasso","avg_elasticnet",
#                  "avg_lasso", "avg_shim","pc_elasticnet",
#                  "pc_lasso", "pc_shim")
# 
# DT.summary[,`:=`(name = factor(name, levels = levels.name, labels = labels.name))]
# 
# DT.summary$name %>% table


# used for boxplots
DT.long2 <- DT.long %>%
  tidyr::unite(name, summary, model)

DT.long2[, table(name)]
levels.name <- c("na_MARS", "avg_MARS", "pc_MARS")
labels.name <- c("MARS","avg_MARS", "pc_MARS")

DT.long2[,`:=`(name = factor(name, levels = levels.name, labels = labels.name))]
DT.long2[, table(name)]
DT.long2[, table(method)]
DT.long2[, table(rho)]
DT.long2[, table(p)]
DT.long2[, table(n)]
DT.long2[, table(nActive)]
DT.long2[, table(cluster_distance)]
DT.long2[, table(Ecluster_distance)]
DT.long2[, table(SNR)]
DT.long2[, table(betaMean)]
DT.long2[, table(alphaMean)]
DT.long2[, table(measure)]
DT.long2[measure=="RMSE"][,boxplot(value)]


DT.long2[cluster_distance=="TOM"][measure=="RMSE"][rho==0.9][SNR==1]$value %>% boxplot.stats()

## ---- tpr-vs-shat ----
# p <- DT.long[measure %in% c("TPR","Shat")][alphaMean==2][SNR==0.2] %>%
#   tidyr::spread(measure, value) %>%
#   tidyr::unite(name, summary, model)
# p[,`:=`(method=factor(method, levels = levels),name = factor(name, levels = levels.name, labels = labels.name) )]
# 
# p %>%
#   ggplot(aes(x = Shat, y = TPR, color=method)) +
#   geom_point(size=2.5, aes(shape=method)) +
#   ylab("True positive rate")+
#   xlab("number of non-zero estimated coefficients")+
#   facet_grid(alphaMean+SNR~rho) + legend.stuff

#ggsave(paste(Sys.getenv("HOME"),"eclust/simulation/simulation1/plots/TPR_vs_Shat.png", sep = "/"))

## ---- tpr-vs-fpr-cor-vs-tom ----

pd <- position_dodge(width = 1)
group.colors <- c(SEPARATE = "#F8766D", CLUST = "#00BA38" , ECLUST = "#619CFF")
appender1 <- function(string) TeX(paste("$\\rho = $", string))
appender2 <- function(string) TeX(paste("$SNR = $", string))
scaleFUN <- function(x) sprintf("%g", x)

tpr_fpr <- copy(DT.long[measure %in% c("TPR","FPR")]) 

tpr_fpr <- tpr_fpr %>%
  tidyr::spread(measure, value) %>%
  tidyr::unite(name, summary, model)

tpr_fpr[, str(method)]
tpr_fpr[, str(name)]
data.table::setnames(tpr_fpr, "name", "model")

# for (i in c(0.2, 1, 2)) {
#   for (j in c(0.5, 2)) {
#     
#     i = 0.2; j=2
#       ggplot(p[alphaMean==j], aes(x = FPR, y = TPR, color = method)) +
#       geom_jitter(size = 2.5, aes(shape = model)) +
#       facet_grid(SNR + rho ~ cluster_distance) +
#       ylab("True positive rate") +
#       xlab("False positive rate") 
#       
#     
#     ggsave(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim2-sept8/tpr_fpr_snr%g_alpha%g_corr_vs_tom_sim2.png",i,j),
#            width = 13, height = 12)
#   }
# }


p1 <- ggplot(tpr_fpr[SNR==0.2],
             aes(x = FPR, y = TPR, color = method)) +
  geom_point(size = 2.5, aes(shape = model)) +
  ylab("True positive rate") +
  xlab("False positive rate") + 
  facet_grid(rho + SNR ~ cluster_distance,
             scales = "fixed",
             labeller = labeller(rho = as_labeller(appender1,
                                                   default = label_parsed),
                                 SNR = as_labeller(appender2,default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) +
  panel_border()



p2 <- ggplot(tpr_fpr[SNR==1],
             aes(x = FPR, y = TPR, color = method)) +
  geom_point(size = 2.5, aes(shape = model)) +
  ylab("") +
  xlab("False positive rate") + 
  facet_grid(rho + SNR  ~ cluster_distance, 
             scales = "fixed", 
             labeller = labeller(rho = as_labeller(appender1, 
                                                   default = label_parsed),
                                 SNR = as_labeller(appender2,default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) +
  panel_border()

p3 <- ggplot(tpr_fpr[SNR==2],
             aes(x = FPR, y = TPR, color = method)) +
  geom_point(size = 2.5, aes(shape = model)) +
  ylab("") +
  xlab("False positive rate") + 
  facet_grid(rho  + SNR ~ cluster_distance, 
             scales = "fixed", 
             labeller = labeller(rho = as_labeller(appender1, 
                                                   default = label_parsed),
                                 SNR = as_labeller(appender2,default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt")) +
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) +
  panel_border()

dev.off()
prow <- plot_grid(p1 + theme(legend.position="none"),
                  p2 + theme(legend.position="none"),
                  p3 + theme(legend.position="none"),
                  NULL,
                  align = 'hv',
                  #ncol = 3,
                  nrow = 2,
                  labels = c(LETTERS[1:3],""),
                  hjust = -1)
prow
grobs <- ggplotGrob(p1 + theme(legend.position="right"))$grobs
legend_b <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
# p <- plot_grid( prow, legend_b, nrow = 2, rel_heights = c(1, .20))

p <- prow + draw_grob(legend_b, 0.3, -0.25)

save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim3-sept27/%s_sim3.png","tpr_fpr",05),
          p, base_aspect_ratio = 1.2, ncol = 2, nrow = 2, base_height = 5)



## ---- tpr-vs-fpr-tom-only ----

pd <- position_dodge(width = 1)
group.colors <- c(SEPARATE = "#F8766D", CLUST = "#00BA38" , ECLUST = "#619CFF")
appender1 <- function(string) TeX(paste("$\\rho = $", string))
appender2 <- function(string) TeX(paste("$SNR = $", string))
scaleFUN <- function(x) sprintf("%g", x)

tpr_fpr <- copy(DT.long[measure %in% c("TPR","FPR")]) 

tpr_fpr <- tpr_fpr %>%
  tidyr::spread(measure, value) %>%
  tidyr::unite(name, summary, model)

tpr_fpr[, str(method)]
tpr_fpr[, str(name)]
data.table::setnames(tpr_fpr, "name", "model")

tpr_fpr$model %>% unique


tpr_fpr[, `:=`(model = factor(model, levels = c("na_MARS","avg_MARS", "pc_MARS"),
                              labels = c("MARS", "avg_MARS", "pc_MARS")))]

for (cdist in c("Correlation","TOM")) {
  p1 <- ggplot(tpr_fpr[cluster_distance==cdist],
               aes(x = FPR, y = TPR, color = method)) +
    geom_point(size = 2.5, aes(shape = model)) +
    ylab("True positive rate") +
    xlab("False positive rate") + 
    facet_grid(SNR ~ rho, 
               scales = "fixed", 
               labeller = labeller(rho = as_labeller(appender1, 
                                                     default = label_parsed),
                                   SNR = as_labeller(appender2,default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position = "bottom") +
    scale_x_continuous(labels=scaleFUN) + 
    scale_y_continuous(labels=scaleFUN) +
    panel_border()
  
  
  save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim3-sept27/tpr_fpr_%s_sim3.png",cdist),
            p1, base_aspect_ratio = 1.3, nrow = 1, base_width = 11, base_height = 11)
  
}

## ---- mse ----

pd <- position_dodge(width = 1)
group.colors <- c(SEPARATE = "#F8766D", CLUST = "#00BA38" , ECLUST = "#619CFF")
appender1 <- function(string) TeX(paste("$\\rho = $", string))
appender2 <- function(string) TeX(paste("$SNR = $", string))


# legend.stuff.2 <- theme_bw()+
#   theme(legend.position = "top")+
#   theme(axis.text.x  = element_text(angle=90, vjust=0.7, size=15),
#         axis.text.y  = element_text(size=15),
#         axis.title.x = element_text(face="bold", colour="#990000", size=15),
#         axis.title.y = element_text(face="bold", colour="#990000", size=20),
#         title = element_text(size=16),legend.text = element_text(colour="blue", size = 16),
#         strip.text = element_text(size=20))+
#   theme(legend.key.width=unit(1, "inches"))



# for (i in c(0.5,2)) {
#  
#     ggplot(DT.long2[measure=="RMSE"][alphaMean==i], 
#            aes(x = name, y = value, fill = method)) + 
#       geom_boxplot(position = pd) +
#       guides(color = guide_legend(title = "method")) +
#       xlab("model") +
#       ylab("RMSE") +
#       facet_grid(SNR ~ rho, scales="free") +
#       legend.stuff.2
#     
#     ggsave(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim2-sept8/RMSE_alpha%g_sim2.png",i),
#            width = 13, height = 12)
#   
# }

DT.long2[,table(measure)]
DT.long2[,table(method)]
DT.long2[,table(name)]
DT.long2[,table(cluster_distance, Ecluster_distance)]

# p1 <- ggplot(DT.long2[cluster_distance == "TOM"][measure=="RMSE"],
#              aes(x = name, y = value, fill = method)) +
#   geom_boxplot(position = pd) +
#   guides(color = guide_legend(title = "method")) +
#   xlab("") +
#   ylab("RMSE") +
#   facet_grid(SNR ~ rho,
#              scales="free",
#              labeller = labeller(rho = as_labeller(appender1,
#                                                    default = label_parsed),
#                                  SNR = as_labeller(appender2,default = label_parsed))) +
#   scale_fill_manual(values = group.colors) +
#   theme(plot.margin = unit(c(6,0,6,0), "pt"),legend.position="bottom") +
#   panel_border()
# 
# p1
for (cdist in c("Correlation","TOM")) {
  
  ll <- DT.long2[cluster_distance == cdist][measure=="RMSE"][, `:=`(ymin = boxplot.stats(value)$stats[1],
                                                                    lower = boxplot.stats(value)$stats[2],
                                                                    middle = boxplot.stats(value)$stats[3],
                                                                    upper = boxplot.stats(value)$stats[4],
                                                                    ymax = boxplot.stats(value)$stats[5]), 
                                                             by = list(method, SNR, rho)]
  
  
  p1 <- ggplot(ll,
               aes(x = name, lower=lower, upper=upper, middle=middle, ymin=ymin,
                   ymax=ymax, fill = method)) +
    geom_boxplot(position = pd, stat="identity") +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("RMSE") +
    facet_grid(SNR ~ rho,
               scales="free",
               labeller = labeller(rho = as_labeller(appender1,
                                                     default = label_parsed),
                                   SNR = as_labeller(appender2,default = label_parsed))) +
    scale_fill_manual(values = group.colors) +
    theme(plot.margin = unit(c(6,0,6,0), "pt"),legend.position="bottom") +
    panel_border()
  
  
  save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim3-sept27/RMSE_%s_sim3.png",cdist),
            p1, base_aspect_ratio = 1.3, nrow = 1, base_width = 11, base_height = 11)
}


## ---- tpr-fpr ---- 

for (cdist in c("Correlation","TOM")) {
  for (m in c("TPR", "FPR")) {
    
    ylabel <- switch(m, 
                     TPR = "True Positive Rate",
                     FPR = "False Positive Rate")
    
    p1 <- ggplot(DT.long2[cluster_distance == cdist][measure==m],
                 aes(x = name, y = value, fill = method)) +
      geom_boxplot(position = pd) +
      guides(color = guide_legend(title = "method")) +
      xlab("") +
      ylab(ylabel) +
      facet_grid(SNR ~ rho,
                 scales="free",
                 labeller = labeller(rho = as_labeller(appender1,
                                                       default = label_parsed),
                                     SNR = as_labeller(appender2,default = label_parsed))) +
      scale_fill_manual(values = group.colors) +
      theme(plot.margin = unit(c(6,0,6,0), "pt"),legend.position="bottom") +
      panel_border()
    
    p1
    
    save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim3-sept27/%s_sim3.png",m),
              p1, base_aspect_ratio = 1.3, nrow = 1, base_width = 11, base_height = 11)
    
  }
}

## ---- jacc-pearson-spearman ----

DT.long2[cluster_distance == "TOM"][measure %in% c("jacc")][is.na(value)]
DT.long2[,table(alphaMean)]

for (cdist in c("Correlation","TOM")) {
  for (m in c("jacc")) {
    
    ylabel <- switch(m, 
                     jacc = "Jaccard Index",
                     spearman = "Spearman Correlation",
                     pearson = "Pearson Correlation")
    
    p1 <- ggplot(DT.long2[cluster_distance == cdist][measure==m],
                 aes(x = name, y = value, fill = method)) +
      geom_boxplot(position = pd) +
      guides(color = guide_legend(title = "method")) +
      xlab("") +
      ylab(ylabel) +
      facet_grid(SNR ~ rho, 
                 scales="fixed", 
                 labeller = labeller(rho = as_labeller(appender1, 
                                                       default = label_parsed),
                                     SNR = as_labeller(appender2,default = label_parsed))) +
      scale_fill_manual(values = group.colors) + 
      # scale_y_continuous(limits = c(0,1)) + 
      theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
      panel_border()
    
    save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim3-sept27/%s_%s_sim3.png",m, cdist),
              p1, base_aspect_ratio = 1.3, nrow = 1, base_width = 11, base_height = 11)
    
  } 
}

## ---- correct-sparsity ----

# ggplot_build(th)$data[[1]]$fill %>% unique()
# pen color is pinkinsh: "#F8766D"
# clust color is green: "#00BA38" 
# eclust color is blue: "#619CFF"

DT.long2[,table(measure)]
DT.long2[, table(name)]
DT.long2[,table(cluster_distance, Ecluster_distance)]


for (cdist in c("Correlation","TOM")) {
  p1 <- ggplot(DT.long2[cluster_distance == cdist][measure=="CorrectSparsity"][name %ni% c("MARS")],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("Correct Sparsity") +
    facet_grid(SNR ~ rho, 
               scales="fixed", 
               labeller = labeller(rho = as_labeller(appender1, 
                                                     default = label_parsed),
                                   SNR = as_labeller(appender2,default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim3-sept27/CorrectSparsity_%s_sim3.png",cdist),
            p1, base_aspect_ratio = 1.3, nrow = 1, base_width = 11, base_height = 11)
  
}
