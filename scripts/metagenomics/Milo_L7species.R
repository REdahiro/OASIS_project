#-------------------------------------------------------------------------------------------------#
# Differential abundance analysis of metagenome using Milo
# dir: /work22/home/redahiro/analysis/OASIS/metagenome
# env: Seurat4.2
# date: 2023.08.10
#-------------------------------------------------------------------------------------------------#

# load library
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(patchwork)
library(tidyverse)

# args
args <- commandArgs(trailingOnly = TRUE)
print(args)
cov = args[1]                         # Each L7 species abundance (Ex: Ruminococcus gnavus)
proportion = as.numeric(args[2])

##--- Load Data
sce <- readRDS("SCE/SCE_HC131_dims_from_ALLPBMC.rds")   # scRNA-seq data of healthy subjects with metagenomics data.

meta <- colData(sce) %>% as_tibble %>% select(cov) %>% setNames("covariate")
sce$covariate <- meta$covariate

print("sce add colData")
sce %>% colData

##--------------------------------------##
##    Differential abundance testing    ##
##--------------------------------------##

##--- Step1: Create a Milo object
milo <- Milo(sce)

##--- Step2: Construct KNN graph
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "PCA")   # K: deault 20 (FindNeighbors)  using the same value for k used for KNN graph building for clustering and UMAP visualization.

##--- Step3: Defining representative neighbourhoods on the KNN graph
milo <- makeNhoods(milo, prop = proportion, k = 30, d = 30, refined = TRUE, reduced_dims = "PCA") 

# check average neighbourhood size: over 5 x N_samples is required 
p.NhSize = plotNhoodSizeHist(milo)
ggsave(file = paste0("./Milo/Distribution_Neighbourhood.size_",cov,".png"), plot = p.NhSize, width = 10, height = 7, dpi = 100)

##--- Step4: Counting cells in neighbourhoods
# adding n(neighbourhoods) * m(experimental samples) matrix
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="SG_ID")
head(nhoodCounts(milo))      # easily can check whether there is unbalanced distribution or not   

##--- Step5: Defining experimental design
design <- data.frame(colData(milo))[,c("SG_ID", "covariate", "Age", "Sex", "Dataset2", "Dataset3")]
design <- distinct(design)
rownames(design) <- design$SG_ID

##--- Step6: Computing neighbourhood connectivity---##
milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")

##----- Step7: Testing -----##
da_results <- testNhoods(milo, design = ~ Age+Sex+Dataset2+Dataset3+covariate, design.df = design)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 

milo <- buildNhoodGraph(milo)
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=0.1)     # alpha: FDR

ggsave(file = paste0("./Milo/Milo_",cov,"_FDR_0.1.png"), plot = nh_graph_pl, width = 8.5, height = 7.5, dpi = 300)
# ggsave(file = "./Milo/Milo_UMPA_Case_vs_HC.pdf", plot = nh_graph_pl, width = 8.5, height = 7.5, dpi = 50)

saveRDS(milo, paste0("./Milo/Milo_",cov,".rds"))

#----------------------------------------------#
#              No FDR filtering                #
#----------------------------------------------# 

##--- Step5: Defining experimental design
design <- data.frame(colData(milo))[,c("SG_ID", "covariate", "Age", "Sex", "Dataset2", "Dataset3")]
design <- distinct(design)
rownames(design) <- design$SG_ID

##--- Step6: Computing neighbourhood connectivity---##
milo <- calcNhoodDistance(milo, d=30, reduced.dim = "PCA")

##----- Step7: Testing -----##
da_results <- testNhoods(milo, design = ~ Age+Sex+Dataset2+Dataset3+covariate, design.df = design)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 
milo <- buildNhoodGraph(milo)
nh_graph_pl <- plotNhoodGraphDA(milo, da_results, layout="UMAP",alpha=1)     # alpha: FDR

ggsave(file = paste0("./Milo/Milo_",cov,"_noFDR.png"), plot = nh_graph_pl, width = 8.5, height = 7.5, dpi = 300)

##-----------------------------------------##
##         downstream analysis             ##
##-----------------------------------------##

#-- Assign a cell type label to each neighbourhood
da_results <- annotateNhoods(milo, da_results, coldata_col = "l3") 

# making order
l3_cluster <- c("CD4_Naive","CD4_TCM","CD4_TEM","CD4_CTL","Treg",
                "CD8_Naive","CD8_TCM","CD8_TEM","CD8_CTL","MAIT","gdT","Pro_T",
                "NK_CD56dim","NK_cytokine","NK_HLA","NK_CD56bright",
                "B_Naive1","B_Naive2","B_Intermediate","B_Memory","B_Activated","PB",
                "cMono_S100A","cMono_IL1B","intMono","ncMono","cDC","pDC")
da_results$l3 <- factor(da_results$l3, levels=rev(l3_cluster))

# Export DA results
write.csv(da_results,paste0("./Milo/DA_results_",i,".csv"),row.names=F)

# https://github.com/MarioniLab/miloR/issues/181
p.DA_celltype = plotDAbeeswarm(da_results, group.by = "l3", alpha=1)   # if png dont work, setting alpha
ggsave(file = paste0("./Milo/BeeswarmPlot_DA_",i,"_noFDR_withMix.pdf"), plot = p.DA_celltype, width = 6, height = 11, dpi = 300)

