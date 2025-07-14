#---------------------------------------------------------------------------------------------------#
# Permutation of edgeR with PRS @l1 levels             
# env: Seurat4.2 
# dir: /work22/home/redahiro/analysis/OASIS/PRS
# date: 2023.12.15
#---------------------------------------------------------------------------------------------------#

library(edgeR)
library(tidyverse)
library(optparse)

option_list = list(
  make_option("--phenotype", type="character", default=FALSE,
              help="Phenotype of GWAS."),
  make_option("--cell_type", type="character", default=FALSE,
              help="Setting target cluster: Mono or DC etc."),
  make_option("--gene_filter", type="integer", default=10,
              help="% expression."),
  make_option("--num_permu", type="integer", default=100,
              help="Number of permutation times."))
opt = parse_args(OptionParser(option_list=option_list))

# Report
print(paste0('cell type: ', opt$cell_type))
cat('\n\n')

phenotype <- opt$phenotype
cell_type <- opt$cell_type
gene_filter <- opt$gene_filter
num_permu <- opt$num_permu
path="/work22/home/redahiro/analysis/OASIS/sceQTL/PseudoBulk/"


#---------------------------------------------------#
#                   Data import                     #
#---------------------------------------------------#

# Phenotype
pheno_base <- read_csv("/work22/home/redahiro/analysis/OASIS/phenotype/Phenotype_COV.info_OASIS_eQTL_cov88_hc146.csv") %>% 
           mutate(ID=SG_ID) %>% select(-SG_ID) %>% select(ID,everything())

# Pseudobulk Count data
count_data <- read_delim(paste0("DEG/Count/PB_perSample_counts_",gene_filter,"_",cell_type,".txt"), "\t", col_names=T)

# sample pick up: cell count>=10 
QCed.sample <- read_csv(paste0(path,cell_type,"/l1/CellCount/",cell_type,".csv")) %>% setNames(c("ID","Count")) %>%
                        filter(Count>=10) %>%　.$ID 
     
# filtering count data 
count_data <- count_data %>% select(Gene,QCed.sample) %>% select(-KF_sample)
analysis_ID <- colnames(count_data) %>% as_tibble %>% setNames("ID")

# PRS: quantile calculated after removing KS sample
prs <- read_delim(paste0("PRScsx/",phenotype,"_rmHLA_META_phi1e-04.profile_edit")," ", col_names=T) %>%
           mutate(ID=IID,PRS=SCORESUM) %>% select(ID,PRS) %>% filter(ID %in% analysis_ID$ID)   
quantile <- prs$PRS %>% quantile()

pheno <- analysis_ID %>% inner_join(prs, by="ID") %>% inner_join(pheno_base, by="ID") %>%
              mutate(PRS_4 = if_else(PRS>quantile[4],4,
                             if_else(PRS>quantile[3],3,
                             if_else(PRS>quantile[2],2,1))))

#-----------------------------------#
#           Permutation             #
#-----------------------------------#
Sig <- list() 

for(per in c(1:num_permu)){
print(paste0("START: permutation ", per))
cate4 <- pheno$PRS_4
set.seed(per)      # set.seed for loop
cate4_p <- sample(cate4)
pheno_p <- pheno %>% mutate(PRS_4 = cate4_p)

#---------------------------------------------------#
#               COVID19 analysis                    #
#---------------------------------------------------#

# pheno_data editting
COV_ID <- pheno_p %>% filter(Status==1) %>% .$ID
COV <- count_data %>% select(Gene,COV_ID) %>% as.data.frame %>% column_to_rownames("Gene")

COV_ID <- colnames(COV) %>% as_tibble %>% setNames("ID")
pheno_cov <- COV_ID %>% inner_join(pheno_p, by="ID")

PRS_4 <- pheno_cov$PRS_4    
Severity <- pheno_cov$Severity 
Age <- pheno_cov$Age
Sex <- pheno_cov$Sex
Version <- pheno_cov$Version

    #----- 4 bin -----# 
    dge <- DGEList(COV, group = PRS_4)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~PRS_4+Severity+Age+Sex+Version)                  
     
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit,coef="PRS_4") 

    tt <- topTags(qlf, n = Inf)$table %>% rownames_to_column("Gene") %>% as_tibble
    n_Sig_COV1 <- tt %>% filter(FDR<0.05) %>% .$Gene %>% length
    n_Sig_COV2 <- tt %>% filter(FDR<0.1) %>% .$Gene %>% length

#---------------------------------------------------#
#               Healthy analysis                    #
#---------------------------------------------------#

# pheno_data editting
HC_ID <- pheno_p %>% filter(Status==0) %>% .$ID
HC <- count_data %>% select(Gene,HC_ID) %>% as.data.frame %>% column_to_rownames("Gene")

HC_ID <- colnames(HC) %>% as_tibble %>% setNames("ID")
pheno_hc <- HC_ID %>% inner_join(pheno_p, by="ID")

PRS_4 <- pheno_hc$PRS_4    
Age <- pheno_hc$Age
Sex <- pheno_hc$Sex
Version <- pheno_hc$Version

    #----- 4 bin -----# 
    dge <- DGEList(HC, group = PRS_4) 
    dge <- calcNormFactors(dge)
    design <- model.matrix(~PRS_4+Age+Sex+Version)                  
     
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit,coef="PRS_4") 

    tt <- topTags(qlf, n = Inf)$table %>% rownames_to_column("Gene") %>% as_tibble
    n_Sig_HC1 <- tt %>% filter(FDR<0.05) %>% .$Gene %>% length
    n_Sig_HC2 <- tt %>% filter(FDR<0.1) %>% .$Gene %>% length   

Sig[[per]] <- tibble(Permu=per,N_Sig_COV_FDR0.05=n_Sig_COV1,N_Sig_COV_FDR0.1=n_Sig_COV2,
                     N_Sig_HC_FDR0.05=n_Sig_HC1,N_Sig_HC_FDR0.1=n_Sig_HC2)

}

Sig_result <- do.call(rbind,Sig)
write.csv(Sig_result, paste0("DEG/",phenotype,"/Permutation_",num_permu,"_edgeR_",cell_type,"_",gene_filter,"_PRS_4bin.csv"), row.names=F)

