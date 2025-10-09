rm(list = ls())
library(dplyr)
library(multiDEGGs)

metadata <- readRDS("../Data/R4RA/metadata.R4RA.RDS")
syn.RNAseq <- readRDS("../Data/R4RA/txi.R4RA.vst.gene.symbols.RDS")
Olink <- readRDS("../Data/R4RA/Olink.raw.R4RA_matched_syn.meta.RDS")
proteomics <- readRDS("../Data/R4RA/protdata.R4RA.RDS")
phosphoproteomics <- readRDS("../Data/R4RA/phosdata.R4RA.RDS")

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
  dplyr::filter(Visit == 3) %>%    # baseline data only 
  dplyr::filter(Randomised.Medication == "Tocilizumab")

deggs_object.multi <- get_diffNetworks(assayData = assayData_list,
                                       metadata = my.metadata,
                                       category_variable = "DAS28.CRP.rem.V9",
                                       percentile_vector = seq(0.25, 0.98, by = 0.05),
                                       regression_method = "rlm",
                                       padj_method = "q.value",
                                       verbose = TRUE,
                                       show_progressBar = T,
                                       cores = 4)

View_diffNetworks(deggs_object.multi)

