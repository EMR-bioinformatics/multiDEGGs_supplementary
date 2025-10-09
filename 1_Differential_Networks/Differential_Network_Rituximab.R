rm(list = ls())
library(dplyr)
library(multiDEGGs)

metadata <- readRDS("../Data/STRAP/STRAP.metadata_shiny041023.RDS")
syn.RNAseq <- readRDS("../Data/STRAP/STRAP.vst.RDS")
Olink <- readRDS("../Data/STRAP/STRAP.olink.RDS")
proteomics <- readRDS("../Data/STRAP/STRAP.mass.spec.proteomics.RDS")
phosphoproteomics <- readRDS("../Data/STRAP/STRAP.mass.spec.phosphoproteomics.RDS")

# Phospho only (for each gene I'll take the average of its phosphorilation sites)
gene_symbols <- gsub("\\([^)]*\\)", "", rownames(phosphoproteomics))
phosphoproteomics <- aggregate(phosphoproteomics, by=list(GeneSymbol=gene_symbols),
                               FUN=function(x) mean(x, na.rm=TRUE))
rownames(phosphoproteomics) <- phosphoproteomics$GeneSymbol
phosphoproteomics$GeneSymbol <- NULL

assayData_list <- list ("syn.RNAseq" = syn.RNAseq, 
                        "proteomics" = proteomics,
                        "phosphoproteomics" = phosphoproteomics,
                        "Olink" = Olink)

my.metadata <- metadata %>%
  dplyr::filter(Visit == 3) %>%  # baseline data only 
  dplyr::filter(Randomised.Medication == "Rituximab")

deggs_object.multi <- get_diffNetworks(assayData = assayData_list,
                                       metadata = my.metadata,
                                       category_variable = "DAS28.CRP.EULARresp.bin.V7",
                                       percentile_vector = seq(0.25, 0.98, by = 0.05),
                                       regression_method = "rlm",
                                       padj_method = "q.value",
                                       verbose = TRUE,
                                       show_progressBar = T,
                                       cores = 4)

View_diffNetworks(deggs_object.multi)

plot_regressions(deggs_object.multi, gene_A = "ITGB3", gene_B = "CRKL", 
                 title = "RNA-seq", legend_position = "bottomright")



