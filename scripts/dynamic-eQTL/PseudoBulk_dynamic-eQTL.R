#---------------------------------------------------------------------------------------------------------------------------#
# Script to run respone-eQTL analysis with Pseudo-bulk approach 
# 
# Ref1: https://github.com/powellgenomicslab/onek1k_phase1/blob/main/lineage_dynamic_analysis/dyn-eqtl_mapping.R (Science 2022)
# Ref2: https://github.com/immunogenomics/sceQTL/tree/main/scripts/singlecell (Nature 2022)
#
# dir: /work22/home/redahiro/analysis/OASIS/sceQTL/reQTL
# env: Seurat4.2
# date: 2024.02.12
#---------------------------------------------------------------------------------------------------------------------------#

# Load libraries
library(lme4)
library(Matrix)
library(data.table)
library(tidyverse)
library(optparse)

options(stringsAsFactors=F)

# options parser
option_list = list(
  make_option("--cluster", type="character", default=FALSE,
              help="Major Six cell types: CD4T or CD8T or NK or B or MD."),
  make_option("--module", type="character", default=FALSE,
              help="Module list of analysis."),
  make_option("--num_cate", type="integer", default=5,
              help="Number of categoies of module scores."))
opt = parse_args(OptionParser(option_list=option_list))

cluster <- opt$cluster
module <- opt$module                 
num_cate <- opt$num_cate

path = "/work22/home/redahiro/analysis/OASIS/sceQTL/reQTL/PseudoBulk/"

##------------------------------##
##        Import data           ##
##------------------------------##

###--- eGene & eSNP paris ---#
# autosomal
eGene_SNP <- read_delim(paste0(path,"tensorQTL/",cluster,"/",module,"/eSNP/",module,"_Cate.",num_cate,"_eSNP_results.txt"),"\t",col_names=T) %>%
               mutate(variant_id_2 = variant_id) %>%
               separate(variant_id_2,c("CHR","REST"), extra="drop") %>%
               filter(!CHR == "chrX") %>% select(gene, variant_id)
eGene_list <- unique(eGene_SNP$gene)

# chrX
eGene_SNP_chrX <- read_delim(paste0(path,"tensorQTL/",cluster,"/",module,"/eSNP/",module,"_Cate.",num_cate,"_eSNP_results.txt"),"\t",col_names=T) %>%
               mutate(variant_id_2 = variant_id) %>%
               separate(variant_id_2,c("CHR","REST"), extra="drop") %>%
               filter(CHR == "chrX") %>% select(gene, variant_id)
eGene_chrX_list <- unique(eGene_SNP_chrX$gene)

###--- Covariates ---###
# Base phenotype
pheno_base <- read_csv("/work22/home/redahiro/analysis/OASIS/phenotype/Phenotype_COV.info_OASIS_eQTL_cov88_hc146.csv") %>%
                      select(SG_ID:Version)

# PC DNA
PC_g <- read_csv("/work22/home/redahiro/analysis/OASIS/phenotype/PCA_genotype_ASA_ALL_cov88_hc148.csv") %>% select(SG_ID,PC1_g,PC2_g)

# PC RNA (Sample * Categories)
PC_r <- read_delim(paste0(path,"Input/",cluster,"/",module,"_Cate.",num_cate,"_Gene.1_Cell.10_qtlInput.Pcs.txt"), "\t", col_names=T) %>%
            mutate(Cate_ID = ...1,
            	     Cate = str_sub(Cate_ID, start=1,end=-16),
                   SG_ID = str_sub(Cate_ID, start=-14,end=-1)) %>%
            select(-Cate_ID) %>%
            select(SG_ID,Cate,PC1:PC15)
# Merge
pheno <- pheno_base %>% left_join(PC_g, by="SG_ID") %>% left_join(PC_r, by= "SG_ID") %>% 
           mutate_at(vars(Cate),as.numeric) %>% mutate(Cate2 = Cate*Cate) %>% select(SG_ID,Cate,Cate2,everything())

###--- Pseudo-bulk expression matrix (Sample * Categories) ---###

# autosomal
exp <- read_delim(paste0(path,"Input/",cluster,"/",module,"_Cate.",num_cate,"_Gene.1_Cell.10_qtlInput.txt"), "\t", col_names=T) %>%
            as.data.frame %>% column_to_rownames("...1") %>% t() %>% as.data.frame %>% rownames_to_column("Cate_ID") %>% as_tibble %>%
            mutate(Cate = str_sub(Cate_ID, start=1,end=-16),
                   SG_ID = str_sub(Cate_ID, start=-14,end=-1)) %>%
            select(-Cate_ID) %>%
            select(SG_ID,Cate,eGene_list) %>% mutate_at(vars(Cate),as.numeric)

# chrX
exp_chrX <- read_delim(paste0(path,"Input/",cluster,"/",module,"_Cate.",num_cate,"_Gene.1_Cell.10_qtlInput.txt"), "\t", col_names=T) %>%
            as.data.frame %>% column_to_rownames("...1") %>% t() %>% as.data.frame %>% rownames_to_column("Cate_ID") %>% as_tibble %>%
            mutate(Cate = str_sub(Cate_ID, start=1,end=-16),
                   SG_ID = str_sub(Cate_ID, start=-14,end=-1)) %>%
            select(-Cate_ID) %>%
            select(SG_ID,Cate,eGene_chrX_list) %>% mutate_at(vars(Cate),as.numeric)

###--- load vcf files ---###

# autosomal
geno <- fread(paste0(path,"tensorQTL/",cluster,"/",module,"/eSNP/vcf/Auto_Cate",num_cate,"_modified.vcf.gz"), header=TRUE, sep='\t', skip='#CHROM') %>%   # data.frame
          select(-c('#CHROM':POS,REF:FORMAT))
colnames_geno <- colnames(geno)
colnames(geno) <- c("variant_id",colnames_geno[2:length(colnames_geno)])

geno <- geno %>% column_to_rownames("variant_id") %>% t() %>% 
            as.data.frame() %>% rownames_to_column("SG_ID") %>% as_tibble 
#geno <- tibble(geno)

# chrX
geno_chrX <- fread(paste0(path,"tensorQTL/",cluster,"/",module,"/eSNP/vcf/chrX_Cate",num_cate,".vcf.gz"), header=TRUE, sep='\t', skip='#CHROM') %>%   # data.frame
          select(-c('#CHROM':POS,REF:FORMAT))
colnames_geno_chrX <- colnames(geno_chrX)
colnames(geno_chrX) <- c("variant_id",colnames_geno_chrX[2:length(colnames_geno_chrX)])

geno_chrX <- geno_chrX %>% column_to_rownames("variant_id") %>% t() %>% 
                 as.data.frame() %>% rownames_to_column("SG_ID") %>% as_tibble 

#####--------------------------------------------------#####
#####            dynamic-eQTL analysis                 #####
#####--------------------------------------------------#####

###  autosomal
out_summary <- list()
out_summary.q <- list()

for(i in c(1:dim(eGene_SNP)[1])){
  
    # pickup gene*variant paris 
    GENE = eGene_SNP[i,1] %>% .$gene
    SNP = eGene_SNP[i,2] %>% .$variant_id

    # creating analysis data.frame
    sub_exp <- exp %>% select(SG_ID,Cate,GENE) %>% setnames(c("SG_ID","Cate","E"))
    sub_geno <- geno %>% select(SG_ID,SNP) %>% setnames(c("SG_ID","G"))

    data <- sub_exp %>% 
              left_join(sub_geno, by="SG_ID") %>%
              left_join(pheno, by=c("SG_ID","Cate")) %>% 
              select(SG_ID,Cate,Cate2,everything()) %>% as.data.frame() 

    data$G <- as.numeric(data$G)
    data$Cate <- as.numeric(data$Cate)
    data$Cate2 <- as.numeric(data$Cate2)
    
    ###--- linear mixed model ---###
    # Fit null model
    fit0 <- lmer(formula=paste0("E ~ G + Cate + Status + Age + Sex + Version + PC1_g + PC2_g +", paste0("PC", 1:15, collapse = " + "), 
                          " + (1 | SG_ID)"), data = data, REML = FALSE)

    # Fit augmented model
    fit1 <- lmer(formula=paste0("E ~ G + Cate + Status + Age + Sex + Version + PC1_g + PC2_g +", paste0("PC", 1:15, collapse = " + "), 
                          " + (1 | SG_ID) + G*Cate"), data = data, REML = FALSE)
    
    model_lrt <- anova(fit0, fit1)

    coef_fit0 <- summary(fit0)$coefficients
    coef_fit1 <- summary(fit1)$coefficients
    out <- data.frame(gene=GENE,snp=SNP,lrt_pval=model_lrt$`Pr(>Chisq)`[2],
                             B_g=coef_fit0[2,1],SE_g=coef_fit0[2,2], T_g=coef_fit0[2,3],
                             B_c=coef_fit0[3,1],SE_c=coef_fit0[3,2], T_c=coef_fit0[3,3],
                             B_gc=coef_fit1[25,1],SE_gc=coef_fit1[25,2], T_gc=coef_fit1[25,3])   
    out_summary[[paste0(GENE,"-",SNP)]] <- out

    ###--- qurdratic mixed model ---###
    # Fit null model
    fit.q0 <- lmer(formula=paste0("E ~ G + Cate + Cate2 + Status + Age + Sex + Version + PC1_g + PC2_g +", paste0("PC", 1:15, collapse = " + "), 
                          " + (1 | SG_ID)"), data = data, REML = FALSE)

    # Fit augmented model
    fit.q1 <- lmer(formula=paste0("E ~ G + Cate + Cate2 + Status + Age + Sex + Version + PC1_g + PC2_g +", paste0("PC", 1:15, collapse = " + "), 
                          " + (1 | SG_ID) + G*Cate + G*Cate2"), data = data, REML = FALSE)
    
    model_lrt.q <- anova(fit.q0, fit.q1)
    out.q <- data.frame(gene=GENE,snp=SNP,lrt_pval=model_lrt.q$`Pr(>Chisq)`[2])
    
    out_summary.q[[paste0(GENE,"-",SNP)]] <- out.q
  }

# summarize list
summary_auto <- do.call(rbind,out_summary) %>% as_tibble %>% arrange(lrt_pval)
summary.q_auto <- do.call(rbind,out_summary.q) %>% as_tibble %>% arrange(lrt_pval)

###---  chrX ---###
rm(out_summary)
rm(out_summary.q)
out_summary <- list()
out_summary.q <- list()

for(i in c(1:dim(eGene_SNP_chrX)[1])){
  
    # pickup gene*variant paris 
    GENE = eGene_SNP_chrX[i,1] %>% .$gene
    SNP = eGene_SNP_chrX[i,2] %>% .$variant_id

    # creating analysis data.frame
    sub_exp <- exp_chrX %>% select(SG_ID,Cate,GENE) %>% setnames(c("SG_ID","Cate","E"))
    sub_geno <- geno_chrX %>% select(SG_ID,SNP) %>% setnames(c("SG_ID","G_base")) %>%
                    mutate(G = if_else(G_base %in% c("0","0|0"), "0",
                               if_else(G_base %in% c("1","1|1"), "2", "1"))) %>%
                    select(SG_ID,G)

    data <- sub_exp %>% 
              left_join(sub_geno, by="SG_ID") %>%
              left_join(pheno, by=c("SG_ID","Cate")) %>%
              as.data.frame() 

    data$G <- as.numeric(data$G)
    data$Cate <- as.numeric(data$Cate)
    data$Cate2 <- as.numeric(data$Cate2)

    ###--- linear mixed model ---###
    # Fit null model
    fit0 <- lmer(formula=paste0("E ~ G + Cate + Status + Age + Sex + Version + PC1_g + PC2_g +", paste0("PC", 1:15, collapse = " + "), 
                          " + (1 | SG_ID)"), data = data, REML = FALSE)

    # Fit augmented model
    fit1 <- lmer(formula=paste0("E ~ G + Cate + Status + Age + Sex + Version + PC1_g + PC2_g +", paste0("PC", 1:15, collapse = " + "), 
                          " + (1 | SG_ID) + G*Cate"), data = data, REML = FALSE)
    
    model_lrt <- anova(fit0, fit1)

    coef_fit0 <- summary(fit0)$coefficients
    coef_fit1 <- summary(fit1)$coefficients
    out <- data.frame(gene=GENE,snp=SNP,lrt_pval=model_lrt$`Pr(>Chisq)`[2],
                             B_g=coef_fit0[2,1],SE_g=coef_fit0[2,2], T_g=coef_fit0[2,3],
                             B_c=coef_fit0[3,1],SE_c=coef_fit0[3,2], T_c=coef_fit0[3,3],
                             B_gc=coef_fit1[25,1],SE_gc=coef_fit1[25,2], T_gc=coef_fit1[25,3])
    
    out_summary[[paste0(GENE,"-",SNP)]] <- out

    ###--- qurdratic mixed model ---###
    # Fit null model
    fit.q0 <- lmer(formula=paste0("E ~ G + Cate + Cate2 + Status + Age + Sex + Version + PC1_g + PC2_g +", paste0("PC", 1:15, collapse = " + "), 
                          " + (1 | SG_ID)"), data = data, REML = FALSE)

    # Fit augmented model
    fit.q1 <- lmer(formula=paste0("E ~ G + Cate + Cate2 + Status + Age + Sex + Version + PC1_g + PC2_g +", paste0("PC", 1:15, collapse = " + "), 
                          " + (1 | SG_ID) + G*Cate + G*Cate2"), data = data, REML = FALSE)
    
    model_lrt.q <- anova(fit.q0, fit.q1)
    out.q <- data.frame(gene=GENE,snp=SNP,lrt_pval=model_lrt.q$`Pr(>Chisq)`[2])
    
    out_summary.q[[paste0(GENE,"-",SNP)]] <- out.q
  }

# summarize list
summary_chrX <- do.call(rbind,out_summary) %>% as_tibble %>% arrange(lrt_pval)
summary.q_chrX <- do.call(rbind,out_summary.q) %>% as_tibble %>% arrange(lrt_pval)

### Making outputdata 
summary_linear <- rbind(summary_auto, summary_chrX) %>%
                      distinct(gene, .keep_all=TRUE) %>%
                      mutate(FDR=p.adjust(lrt_pval, method= "BH")) %>%
                      select(gene,snp,FDR,lrt_pval,everything())

write.table(summary_linear, paste0(path,"tensorQTL/",cluster,"/",module,"/Results/summary_linear_",module,"_Cate.",num_cate,"_update.txt"), quote=F, sep="\t", row.names=F, col.names=T)

summary_quadratic <- rbind(summary.q_auto, summary.q_chrX) %>%
                      distinct(gene, .keep_all=TRUE) %>%
                      mutate(FDR=p.adjust(lrt_pval, method= "BH")) %>%
                      select(gene,snp,FDR,lrt_pval,everything())

write.table(summary_quadratic, paste0(path,"tensorQTL/",cluster,"/",module,"/Results/summary_quadratic_",module,"_Cate.",num_cate,"_update.txt"), quote=F, sep="\t", row.names=F, col.names=T)
