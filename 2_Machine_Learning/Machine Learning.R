# Fit models ---------------------------------------------------------------
rm(list = ls())
library(caret)
library(glmnet)
library(gbm)
library(xgboost)
library(ggplot2)
library(dplyr)
library(randomForest)
library(mda)
library(pls)
library(RhpcBLASctl)
library(devtools)
library(parallel)
library(dplyr)
library(nestedcv)
library(multiDEGGs)
source("compare_mod.R")

n_filter = 40  # 40 for Tocilizumab, 50 for Rituximab (larger sample size)
method = "lm"
sig = "qvalue"
gene.set = "coding" # filter for protein coding genes 
dataset = "R4RA"  # R4RA, STRAP
drug = "Tocilizumab"   # Tocilizumab, Rituximab
outcome.var = "DAS28.CRP.rem.V9"  # DAS28.CRP.rem.V9 (R4RA Toc), DAS28.CRP.EULARresp.bin.V7 (STRAP Rtx)

mods <- list(
  glmnet = list(model="glmnet", args=list(family="binomial",
                                          alphaSet = seq(0.7, 1, 0.1),
                                          min_1se = 0)),
  pls = list(model = "pls"),
  pda = list(model="pda"),
  svmPoly = list(model="svmPoly"),
  gbm = list(model="gbm"),
  rf = list(model="rf"),
  xgbTree = list(model="xgbTree"),
  xgbLinear = list(model="xgbLinear")
)


grd <- list(
  modifyX = "multiDEGGs_filter",
  filterFUN = "ttest_filter",
  filterFUN = "wilcoxon_filter",
  filterFUN = "ranger_filter",
  filterFUN = "glmnet_filter",
  filterFUN = "pls_filter"
)

metadata <- switch (dataset,
                    "R4RA" = readRDS("../Data/R4RA/metadata.R4RA.RDS"),
                    "STRAP" = readRDS("../Data/STRAP/STRAP.metadata_shiny041023.RDS"),
                    "R4RA.STRAP" = readRDS("../Data/R4RA_STRAP/metadata.R4RA.STRAP_ES.RDS")
)

vst <- switch (dataset,
               "R4RA" = readRDS("../Data/R4RA/txi.R4RA.vst.gene.symbols.RDS"),
               "STRAP" = readRDS("../Data/STRAP/STRAP.vst.RDS"),
               "R4RA.STRAP" = readRDS("../Data/R4RA_STRAP/vst.R4RA.STRAP.comBat.RDS")
)

my.metadata <- metadata %>%
  dplyr::filter(Visit == 3) %>%
  dplyr::filter(Randomised.Medication == drug)
my.vst <- as.data.frame(vst[, rownames(my.metadata)])

# filter for protein coding genes 
if (gene.set == "coding"){
  coding.genes <- read.csv("../Data/gene_type.txt", stringsAsFactors = FALSE)[, c(1,2)] %>%
    filter(Gene.type == "protein_coding") %>%
    pull(Gene.name)
  
  coding.genes <- coding.genes[coding.genes %in% rownames(my.vst)]
  my.vst <- my.vst[coding.genes, ]
}

y <- as.factor(my.metadata[, outcome.var])
x <- as.matrix(t(my.vst))
x <- x[,apply(x, 2, var, na.rm=TRUE) != 0] # remove zero variance cols

y <- as.factor(gsub("Non-responder", "Non.responder", y))
y <- as.factor(gsub("Good/mod-responder", "Good.Mod.Responder", y))
y <- as.factor(gsub("Non-remission", "Non.remission", y))

fit <- compare_mod(y, 
                   x,
                   mods,
                   grid = grd,
                   filter_options = list(nfilter = n_filter),
                   modifyX_options = list(
                     keep_single_genes = TRUE,
                     nfilter = n_filter
                   ),
                   modifyX_useY = TRUE,
                   repeats = 16,
                   n_outer_folds = 5,
                   rep.cores = 16
)

saveRDS(fit, paste0("res2/fit.", 
                    paste(dataset, drug, outcome.var, n_filter, method, gene.set,
                          "RDS", sep = ".")))


# Prepare plotting data ---------------------------------------
## comparing multiDEGGs_filter with other statistical filters

dataset = "R4RA"  # R4RA, STRAP
drug = "Tocilizumab"   # Tocilizumab, Rituximab
outcome.var = "DAS28.CRP.rem.V9"  # DAS28.CRP.rem.V9 (R4RA Toc), DAS28.CRP.EULARresp.bin.V7 (STRAP Rtx)
n_filter = 40 # 40 for Tocilizumab, 50 for Rituximab (larger sample size)
method = "lm"
sig = "qvalue"
gene.set = "coding"

fit <- readRDS(paste0("res2/fit.", 
                      paste(dataset, drug, outcome.var, n_filter,
                            method, gene.set,
                            "RDS", sep = ".")))
result <- fit$result
w <- unlist(lapply(result, function(i) is.null(i$res)))
result <- result[w]

st <- do.call(rbind, result)
st$filter <- gsub("none", "multiDEGGs filter", st$filter)
st$filter <- gsub("ttest", "t-test", st$filter)
st$filter <- gsub("ranger", "rf", st$filter)
st$filter <- gsub("_", " ", st$filter)
st$filter <- factor(st$filter, levels = c("multiDEGGs filter", "glmnet filter",
                                          "pls filter", "rf filter",
                                          "t-test filter", "wilcoxon filter"))

saveRDS(st, paste0("res2/st.",
                   paste(dataset, drug, outcome.var, n_filter, method, sig,
                         "RDS", sep = ".")))

# Plot 1 --------------------------------------------------------------------
## comparing multiDEGGs_filter with other statistical filters
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(paletteer)
library(tidyverse)

dataset = "R4RA"  # R4RA, STRAP
drug = "Tocilizumab"   # Tocilizumab, Rituximab
outcome.var = "DAS28.CRP.EULARresp.bin.V7"  # DAS28.CRP.rem.V9 (R4RA Toc), DAS28.CRP.EULARresp.bin.V7 (STRAP Rtx)
n_filter = 50 # 40 for Tocilizumab, 50 for Rituximab (larger sample size)
method = "lm"
sig = "qvalue"
gene.set = "coding"
order = FALSE

st <- readRDS(paste0("res2/st.",
                     paste(dataset, drug, outcome.var, n_filter, method, sig,
                           "RDS", sep = ".")))

# order filters by their decreasing AUC mean value, keep multiDEGGs as first
if(order) {
  st <- st %>%
    mutate(filter = fct_reorder(filter, AUC, .fun = mean, .desc = TRUE)) %>%
    mutate(filter = fct_relevel(filter, "multiDEGGs filter"))
}


p <- ggplot(st, aes(x = filter, y = AUC, col = filter)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.05)) +
  geom_boxplot(outlier.shape = NA, alpha = 0) +
  scale_color_paletteer_d("ggsci::default_nejm", dynamic = F) +
  facet_grid(. ~ model, scales = "free_x", space = "free") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(col = guide_legend(title = NULL)) 
  # labs(title = paste(dataset, drug, n_filter, method, sig))
p

ggsave(plot = p, path = "res2", 
       filename = paste(dataset, drug, outcome.var, n_filter, method, sig, "svg",
                        sep = "."), 
       width = 7.68, height = 4.8, units = "in")


# compute 95% CI 95% for each model-filter conmbination
confidence_intervals <- st %>%
  group_by(model, filter) %>%
  summarise(
    mean_AUC = mean(AUC, na.rm = TRUE),
    median_AUC = median(AUC, na.rm = TRUE),
    CI_lower = quantile(AUC, 0.025, na.rm = TRUE),
    CI_upper = quantile(AUC, 0.975, na.rm = TRUE),
    sd_AUC = sd(AUC, na.rm = TRUE),
    n_repeats = n()
  )


# write csv for supplementary tables 
write.csv(confidence_intervals, 
          paste0("res2/CI.",
                 paste(dataset, drug, outcome.var, n_filter, method, sig, "csv",
                       sep = ".")))

write.csv(st, 
          paste0("res2/st.",
                 paste(dataset, drug, outcome.var, n_filter, method, sig, "csv",
                       sep = ".")))
