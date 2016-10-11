##################################
# R source code file for combining analysis to create plots that contain 
# results from all three simulations
# 
# Created by Sahir, October 1, 2016
# Updated: 
##################################



# Load data ---------------------------------------------------------------

rm(list = ls())
source("~/git_repositories/eclust-simulation-aug2016/packages.R")
source("~/git_repositories/eclust-simulation-aug2016/functions.R")
options(digits = 2, scipen=999)
options(warning.length = 8170)

# Simulation 1 ------------------------------------------------------------

col.names <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/colnames_stab_hydra-sim1-sept10.txt", header = F)$V1

DT <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/sim1-hydra-results-p5000-sept10", stringsAsFactors = FALSE) %>%
  setnames(col.names)

DT[, `:=`(simulation = 1:nrow(DT))]

# this still has all the raw data, but melted

DT.long <- DT %>%
  reshape2::melt(id.vars = c("simulation",colnames(DT)[1:20])) %>%
  tidyr::separate(variable, c("method", "summary", "model", "interaction", "measure"), convert = T) %>%
  as.data.table

DT.long[, `:=`(method = factor(method, levels = c("pen", "clust", "Eclust"), labels = c("SEPARATE", "CLUST", "ECLUST")))]
DT.long[, `:=`(cluster_distance = factor(cluster_distance, levels = c("corr","tom"), labels = c("Correlation", "TOM")))]

# used for boxplots
DT.long2 <- DT.long %>%
  tidyr::unite(name, summary, model)

levels.name <- c("na_elasticnet", "na_lasso","avg_elasticnet",
                 "avg_lasso", "pc_elasticnet",
                 "pc_lasso")

labels.name <- c("enet", "lasso","avg_enet",
                 "avg_lasso", "pc_enet",
                 "pc_lasso")

DT.long2[,`:=`(name = factor(name, levels = levels.name, labels = labels.name))]


DT.long.sim1 <- copy(DT.long)
DT.long2.sim1 <- copy(DT.long2)
rm(DT,DT.long,DT.long2,col.names)


DT.long.sim1[, `:=`(sim = 1)]
DT.long2.sim1[, `:=`(sim = 1)]



# Simulation 2 ------------------------------------------------------------

col.names <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/colnames_stab_hydra-sim2-sept8.txt", header = F)$V1

# DT <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/sim-hydra-results-p3000-no-stability", stringsAsFactors = FALSE) %>%
#   setnames(col.names)

DT <- fread("~/git_repositories/eclust-simulation-aug2016/hydra/results/sim2-hydra-results-p5000-sept8", stringsAsFactors = FALSE) %>%
  setnames(col.names)

DT[, `:=`(simulation = 1:nrow(DT))]

# this still has all the raw data, but melted
options(warning.length = 8170)
DT.long <- DT %>%
  reshape2::melt(id.vars = c("simulation",colnames(DT)[1:20])) %>%
  tidyr::separate(variable, c("method", "summary", "model", "interaction", "measure"), convert = T) %>%
  as.data.table


DT.long[, `:=`(method = factor(method, levels = c("pen", "clust", "Eclust"), labels = c("SEPARATE", "CLUST", "ECLUST")))]
DT.long[, `:=`(cluster_distance = factor(cluster_distance, levels = c("corr","tom"), labels = c("Correlation", "TOM")))]

# used for boxplots
DT.long2 <- DT.long %>%
  tidyr::unite(name, summary, model)

DT.long2[, table(name)]

levels.name <- c("na_elasticnet", "na_lasso","avg_elasticnet",
                 "avg_lasso", "pc_elasticnet",
                 "pc_lasso")

labels.name <- c("enet", "lasso","avg_enet",
                 "avg_lasso", "pc_enet",
                 "pc_lasso")

DT.long2[,`:=`(name = factor(name, levels = levels.name, labels = labels.name))]


DT.long.sim2 <- copy(DT.long)
DT.long2.sim2 <- copy(DT.long2)
rm(DT,DT.long,DT.long2,col.names)


DT.long.sim2[, `:=`(sim = 2)]
DT.long2.sim2[, `:=`(sim = 2)]

DT.long.sim2 <- DT.long.sim2[alphaMean==2]
DT.long2.sim2 <- DT.long2.sim2[alphaMean==2]


# Simulation 3 ------------------------------------------------------------

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

# used for boxplots
DT.long2 <- DT.long %>%
  tidyr::unite(name, summary, model)

DT.long2[, table(name)]
levels.name <- c("na_MARS", "avg_MARS", "pc_MARS")
labels.name <- c("MARS","avg_MARS", "pc_MARS")
DT.long2[,`:=`(name = factor(name, levels = levels.name, labels = labels.name))]


DT.long.sim3 <- copy(DT.long)
DT.long2.sim3 <- copy(DT.long2)
rm(DT,DT.long,DT.long2,col.names)


DT.long.sim3[, `:=`(sim = 3)]
DT.long2.sim3[, `:=`(sim = 3)]

DT.long2.sim3[cluster_distance=="TOM"][measure=="RMSE"][rho==0.9][SNR==1]$value %>% boxplot.stats()

# Combine data ------------------------------------------------------------

all(colnames(DT.long.sim1) %in% colnames(DT.long.sim2) )
all(colnames(DT.long.sim1) %in% colnames(DT.long.sim3) )
all(colnames(DT.long.sim2) %in% colnames(DT.long.sim3) )
all(colnames(DT.long.sim2) %in% colnames(DT.long.sim1) )
all(colnames(DT.long.sim3) %in% colnames(DT.long.sim1) )
all(colnames(DT.long.sim3) %in% colnames(DT.long.sim2) )

all(colnames(DT.long2.sim1) %in% colnames(DT.long2.sim2) )
all(colnames(DT.long2.sim1) %in% colnames(DT.long2.sim3) )
all(colnames(DT.long2.sim2) %in% colnames(DT.long2.sim3) )
all(colnames(DT.long2.sim2) %in% colnames(DT.long2.sim1) )
all(colnames(DT.long2.sim3) %in% colnames(DT.long2.sim1) )
all(colnames(DT.long2.sim3) %in% colnames(DT.long2.sim2) )


DT.long <- rbindlist(list(DT.long.sim1,DT.long.sim2, DT.long.sim3))
DT.long2 <- rbindlist(list(DT.long2.sim1,DT.long2.sim2, DT.long2.sim3))



## ---- tpr-vs-fpr ----

pd <- position_dodge(width = 1)
group.colors <- c(SEPARATE = "#F8766D", CLUST = "#00BA38" , ECLUST = "#619CFF")
appender1 <- function(string) TeX(paste("Simulation ", string))
appender2 <- function(string) TeX(paste("$SNR = $", string))
scaleFUN <- function(x) sprintf("%g", x)

tpr_fpr <- copy(DT.long[measure %in% c("TPR","FPR")]) 

tpr_fpr <- tpr_fpr %>%
  tidyr::spread(measure, value) %>%
  tidyr::unite(name, summary, model)

tpr_fpr[, str(method)]
tpr_fpr[, str(name)]
data.table::setnames(tpr_fpr, "name", "model")

p1 <- ggplot(tpr_fpr[SNR==1][rho==0.9][cluster_distance=="TOM"],
             aes(x = FPR, y = TPR, color = method)) +
  # geom_point(size = 2.5, aes(shape = model)) +
  geom_point(size = 2.5, aes(shape = method)) +
  ylab("True positive rate") +
  xlab("False positive rate") + 
  facet_grid(~sim, 
             scales = "fixed", 
             labeller = labeller(sim = as_labeller(appender1, 
                                                   default = label_parsed),
                                 SNR = as_labeller(appender2,default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  # theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") +
  theme(legend.position="bottom") +
  scale_x_continuous(labels=scaleFUN) + 
  scale_y_continuous(labels=scaleFUN) +
  panel_border()

save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim1-2-3-combined/%s_sim123.png","tpr_fpr",05),
          p1, base_width = 8)



## ---- mse ----

pd <- position_dodge(width = 1)
group.colors <- c(SEPARATE = "#F8766D", CLUST = "#00BA38" , ECLUST = "#619CFF")
appender1 <- function(string) TeX(paste("Simulation ", string))


# DT.long2[SNR==1][rho==0.9][cluster_distance=="TOM"][measure=="RMSE"][sim==3]$value %>% boxplot.stats

p1 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance=="TOM"][measure=="RMSE"],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("RMSE") +
  facet_grid( ~ sim, 
             scales="free", 
             labeller = labeller(sim = as_labeller(appender1, 
                                                   default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"),legend.position="bottom") +
  panel_border()

save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim1-2-3-combined/%s_sim123.png","RMSE"),
          p1, nrow = 1, base_width = 18, base_height = 7)

DT.long2[SNR==2][rho==0.9][cluster_distance=="TOM"][measure=="RMSE"][sim==3]

## ---- correct-sparsity ----

# ggplot_build(th)$data[[1]]$fill %>% unique()
# pen color is pinkinsh: "#F8766D"
# clust color is green: "#00BA38" 
# eclust color is blue: "#619CFF"

DT.long2[,table(measure)]
DT.long2[, table(name)]
DT.long2[,table(cluster_distance, Ecluster_distance)]



p1 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="CorrectSparsity"][name %ni% c("enet", "lasso")][sim==1],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("Correct Sparsity") +
  facet_grid( ~ sim, 
             scales="fixed", 
             labeller = labeller(sim = as_labeller(appender1, 
                                                   default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
  panel_border()

 p2 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="CorrectSparsity"][name %ni% c("enet", "lasso")][sim==2],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("") +
  facet_grid( ~ sim, 
              scales="fixed", 
              labeller = labeller(sim = as_labeller(appender1, 
                                                    default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
  panel_border()


 p3 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="CorrectSparsity"][name %ni% c("MARS")][sim==3],
              aes(x = name, y = value, fill = method)) +
   geom_boxplot(position = pd) +
   guides(color = guide_legend(title = "method")) +
   xlab("") +
   ylab("") +
   facet_grid( ~ sim, 
               scales="fixed", 
               labeller = labeller(sim = as_labeller(appender1, 
                                                     default = label_parsed))) +
   scale_fill_manual(values = group.colors) + 
   scale_y_continuous(limits = c(0,1)) + 
   theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
   panel_border()
 
 
prow <- plot_grid(p1 + theme(legend.position="none"),
                  p2 + theme(legend.position="none"),
                  p3 + theme(legend.position="none"),
                  align = 'hv',
                  nrow = 1,
                  labels = c("", "",""),
                  hjust = -1)

grobs <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs
legend_b <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1,.03))
p

save_plot("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim1-2-3-combined/CorrectSparsity_sim123.png",
          p, base_aspect_ratio = 1.3, nrow = 1, base_width = 18, base_height = 7)


## ---- jaccard ----

# ggplot_build(th)$data[[1]]$fill %>% unique()
# pen color is pinkinsh: "#F8766D"
# clust color is green: "#00BA38" 
# eclust color is blue: "#619CFF"

DT.long2[,table(measure)]
DT.long2[, table(name)]
DT.long2[,table(cluster_distance, Ecluster_distance)]


p1 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="jacc"][sim==1],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("Jaccard Index") +
  facet_grid( ~ sim, 
              scales="fixed", 
              labeller = labeller(sim = as_labeller(appender1, 
                                                    default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
  panel_border()

p2 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="jacc"][sim==2],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("") +
  facet_grid( ~ sim, 
              scales="fixed", 
              labeller = labeller(sim = as_labeller(appender1, 
                                                    default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
  panel_border()


p3 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="jacc"][sim==3],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("") +
  facet_grid( ~ sim, 
              scales="fixed", 
              labeller = labeller(sim = as_labeller(appender1, 
                                                    default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
  panel_border()


prow <- plot_grid(p1 + theme(legend.position="none"),
                  p2 + theme(legend.position="none"),
                  p3 + theme(legend.position="none"),
                  align = 'hv',
                  nrow = 1,
                  labels = c("", "",""),
                  hjust = -1)

grobs <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs
legend_b <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
p <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1,.03))
p

save_plot("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim1-2-3-combined/jacc_sim123.png",
          p, base_aspect_ratio = 1.3, nrow = 1, base_width = 18, base_height = 7)



## ---- pearson-spearman ----

# ggplot_build(th)$data[[1]]$fill %>% unique()
# pen color is pinkinsh: "#F8766D"
# clust color is green: "#00BA38" 
# eclust color is blue: "#619CFF"

DT.long2[,table(measure)]
DT.long2[, table(name)]
DT.long2[,table(cluster_distance, Ecluster_distance)]


for (m in c("pearson","spearman")) {
  
  ylabel <- switch(m, 
                   jacc = "Jaccard Index",
                   spearman = "Spearman Correlation",
                   pearson = "Pearson Correlation")
  
  p1 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure==m][sim==1],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab(ylabel) +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  p2 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure==m][sim==2],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("") +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  
  prow <- plot_grid(p1 + theme(legend.position="none"),
                    p2 + theme(legend.position="none"),
                    align = 'hv',
                    nrow = 1,
                    labels = c("", ""),
                    hjust = -1)
  
  grobs <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs
  legend_b <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1,.03))
  p
  
  save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim1-2-3-combined/%s_sim123.png",m),
            p, base_aspect_ratio = 1.3, nrow = 1, base_width = 18, base_height = 7)
}



## ---- stability combined ----
# ggplot_build(th)$data[[1]]$fill %>% unique()
# pen color is pinkinsh: "#F8766D"
# clust color is green: "#00BA38" 
# eclust color is blue: "#619CFF"

pd <- position_dodge(width = 1)
DT.long2[,table(measure)]
DT.long2[, table(name)]
DT.long2[,table(cluster_distance, Ecluster_distance)]
pd <- position_dodge(width = 1)
group.colors <- c(SEPARATE = "#F8766D", CLUST = "#00BA38" , ECLUST = "#619CFF")
appender1 <- function(string) TeX(paste("Simulation ", string))
appender2 <- function(string) TeX(paste("$SNR = $", string))
scaleFUN <- function(x) sprintf("%g", x)

p1 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="jacc"][sim==1],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("Jaccard Index") +
  facet_grid( ~ sim, 
              scales="fixed", 
              labeller = labeller(sim = as_labeller(appender1, 
                                                    default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
  panel_border()

p2 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="jacc"][sim==2],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("") +
  facet_grid( ~ sim, 
              scales="fixed", 
              labeller = labeller(sim = as_labeller(appender1, 
                                                    default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
  panel_border()


p3 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="jacc"][sim==3],
             aes(x = name, y = value, fill = method)) +
  geom_boxplot(position = pd) +
  guides(color = guide_legend(title = "method")) +
  xlab("") +
  ylab("") +
  facet_grid( ~ sim, 
              scales="fixed", 
              labeller = labeller(sim = as_labeller(appender1, 
                                                    default = label_parsed))) +
  scale_fill_manual(values = group.colors) + 
  scale_y_continuous(limits = c(0,1)) + 
  theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
  panel_border()


p_jacc <- plot_grid(p1 + theme(legend.position="none"),
                  p2 + theme(legend.position="none"),
                  p3 + theme(legend.position="none"),
                  align = 'hv',
                  nrow = 1,
                  labels = c("", "",""),
                  hjust = -1)


m <- "pearson"
  
  ylabel <- switch(m, 
                   jacc = "Jaccard Index",
                   spearman = "Spearman Correlation",
                   pearson = "Pearson Correlation")
  
  p1 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure==m][sim==1],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab(ylabel) +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  p2 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure==m][sim==2],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("") +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  
  p_pearson <- plot_grid(p1 + theme(legend.position="none"),
                    p2 + theme(legend.position="none"),
                    align = 'hv',
                    nrow = 1,
                    labels = c("", ""),
                    hjust = -1)
  
m <- "spearman"
  
  ylabel <- switch(m, 
                   jacc = "Jaccard Index",
                   spearman = "Spearman Correlation",
                   pearson = "Pearson Correlation")
  
  p1 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure==m][sim==1],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab(ylabel) +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  p2 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure==m][sim==2],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("") +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  
  p_spearman <- plot_grid(p1 + theme(legend.position="none"),
                    p2 + theme(legend.position="none"),
                    align = 'hv',
                    nrow = 1,
                    labels = c("", ""),
                    hjust = -1)
  
  grobs <- ggplotGrob(p1 + theme(legend.position="bottom"))$grobs
  legend_b <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(p_jacc, p_pearson, p_spearman, legend_b, nrow = 4, rel_heights = c(1,1,1,.10), labels = c("A","B","C",""))
  p
  
  save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim1-2-3-combined/%s_sim123.png","stability"),
            p, base_aspect_ratio = 1.3, nrow = 1, base_width = 16, base_height = 14)


  ## ---- model-fit-combined ----
  
  pd <- position_dodge(width = 1)
  group.colors <- c(SEPARATE = "#F8766D", CLUST = "#00BA38" , ECLUST = "#619CFF")
  appender1 <- function(string) TeX(paste("Simulation ", string))
  appender2 <- function(string) TeX(paste("$SNR = $", string))
  scaleFUN <- function(x) sprintf("%g", x)
  
  tpr_fpr <- copy(DT.long[measure %in% c("TPR","FPR")]) 
  
  tpr_fpr <- tpr_fpr %>%
    tidyr::spread(measure, value) %>%
    tidyr::unite(name, summary, model)
  
  tpr_fpr[, str(method)]
  tpr_fpr[, str(name)]
  data.table::setnames(tpr_fpr, "name", "model")
  
  p_tpr_fpr <- ggplot(tpr_fpr[SNR==1][rho==0.9][cluster_distance=="TOM"],
               aes(x = FPR, y = TPR, color = method)) +
    # geom_point(size = 2.5, aes(shape = model)) +
    geom_point(size = 2.5, aes(shape = method)) +
    ylab("True positive rate") +
    xlab("False positive rate") + 
    facet_grid(~sim, 
               scales = "fixed", 
               labeller = labeller(sim = as_labeller(appender1, 
                                                     default = label_parsed),
                                   SNR = as_labeller(appender2,default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    # theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") +
    theme(legend.position="none") +
    scale_x_continuous(labels=scaleFUN) + 
    scale_y_continuous(labels=scaleFUN) +
    panel_border()
  
  
  p1 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="CorrectSparsity"][name %ni% c("enet", "lasso")][sim==1],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("Correct Sparsity") +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  p2 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="CorrectSparsity"][name %ni% c("enet", "lasso")][sim==2],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("") +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  
  p3 <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance == "TOM"][measure=="CorrectSparsity"][name %ni% c("MARS")][sim==3],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("") +
    facet_grid( ~ sim, 
                scales="fixed", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    scale_y_continuous(limits = c(0,1)) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"), legend.position="bottom") + 
    panel_border()
  
  
  p_correct_sparsity <- plot_grid(p1 + theme(legend.position="none"),
                    p2 + theme(legend.position="none"),
                    p3 + theme(legend.position="none"),
                    align = 'hv',
                    nrow = 1,
                    labels = c("", "",""),
                    hjust = -1)
  
  pd <- position_dodge(width = 1)
  group.colors <- c(SEPARATE = "#F8766D", CLUST = "#00BA38" , ECLUST = "#619CFF")
  appender1 <- function(string) TeX(paste("Simulation ", string))
  
  p_rmse <- ggplot(DT.long2[SNR==1][rho==0.9][cluster_distance=="TOM"][measure=="RMSE"],
               aes(x = name, y = value, fill = method)) +
    geom_boxplot(position = pd) +
    guides(color = guide_legend(title = "method")) +
    xlab("") +
    ylab("RMSE") +
    facet_grid( ~ sim, 
                scales="free", 
                labeller = labeller(sim = as_labeller(appender1, 
                                                      default = label_parsed))) +
    scale_fill_manual(values = group.colors) + 
    theme(plot.margin = unit(c(6,0,6,0), "pt"),legend.position="bottom") +
    panel_border()
  
  p <- plot_grid(p_tpr_fpr, p_correct_sparsity, p_rmse, nrow = 3, rel_heights = c(1,1,1), labels = c("A","B","C"))
  p
  
  save_plot(sprintf("~/git_repositories/eclust-simulation-aug2016/hydra/results/figures/sim1-2-3-combined/%s_sim123.png","modelfit"),
            p, base_aspect_ratio = 1.3, nrow = 1, base_width = 16, base_height = 14)
  