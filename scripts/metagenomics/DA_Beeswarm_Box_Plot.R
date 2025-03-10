#-------------------------------------------------------------------------------------------------#
# Beamswarm plot with boxplot 
# dir: /work22/home/redahiro/analysis/OASIS/metagenome
# env: RASTR
# date:  2024.02.20 
#-------------------------------------------------------------------------------------------------#

# load library
library(ggrastr)
library(tidyverse)
library(ggbeeswarm)

# setting
i = "RG_i"         # L7_species
alpha = 1
group_by = "l3"

# l3 order
l3_cluster <- c("CD4_Naive","CD4_TCM","CD4_TEM","CD4_CTL","Treg",
                "CD8_Naive","CD8_TCM","CD8_TEM","CD8_CTL","MAIT","gdT","Pro_T",
                "NK_CD56dim","NK_cytokine","NK_HLA","NK_CD56bright",
                "B_Naive1","B_Naive2","B_Intermediate","B_Memory","B_Activated","PB",
                "cMono_IL1B","cMono_S100A","intMono","ncMono","cDC","pDC")

# Import of DA results
da_results <- read_csv(paste0("./Milo/DA_results_",i,".csv"))
da_results$l3 <- factor(da_results$l3, levels=l3_cluster)
n_groups <- unique(da_results$l3) %>% length()

p <- ggplot(da_results, aes(x=l3, y=logFC, color=logFC)) +
       geom_hline(yintercept=0, linetype = 1, size = 0.5, color = "black") +
       geom_quasirandom_rast(width=0.3) +
       scale_color_gradient2() +
       guides(color="none") +
       xlab("") + ylab("Log2 Fold Change") +
       geom_boxplot(alpha=0.1,outlier.shape=NA,width=0.25) +
       #coord_flip() +
       # theme_bw(base_size=22) +
       theme_classic(base_size=18) +
       theme(axis.text.x = element_text(size=14, angle=90, hjust=1, vjust=1))

ggsave(file = paste0("./Milo/Beeswarm_Box_DA_",i,"_noFDR_narrow.pdf"), plot = p, width = 11, height = 5, dpi = 300)


