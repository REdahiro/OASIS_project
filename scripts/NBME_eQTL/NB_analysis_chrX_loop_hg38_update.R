#---------------------------------------------------------------------------------------------------------------------------#
# Script to run cell-state-dependent eQTL analysis on only chrX using NB model at single-cell resolution
# 
# Ref1: https://github.com/powellgenomicslab/onek1k_phase1/blob/main/lineage_dynamic_analysis/dyn-eqtl_mapping.R (Science 2022)
# Ref2: https://github.com/immunogenomics/sceQTL/tree/main/scripts/singlecell (Nature 2022)
#
# dir: /work22/home/redahiro/analysis/OASIS/sceQTL/reQTL
# env: Seurat4.2
# date: 2023.10.08
#---------------------------------------------------------------------------------------------------------------------------#

# Load libraries
library(lme4)
library(Matrix)
library(data.table)
library(tidyverse)
library(optparse)
# library(qvalue)

options(stringsAsFactors=F)

# options parser
option_list = list(
  make_option("--cluster", type="character", default=FALSE,
              help="Major Six cell types: CD4T or CD8T or NK or B or MD."),
  make_option("--module", type="character", default=FALSE,
              help="Module list of analysis."),
  make_option("--num_cate", type="integer", default=5,
              help="Number of categoies of module scores."),
  make_option("--test", type="character", default=FALSE,
              help="Pseudobulk dynamic-eQTL analysis: linear or quadratic."))
opt = parse_args(OptionParser(option_list=option_list))

cluster <- opt$cluster
module <- opt$module                 
num_cate <- opt$num_cate
test <- opt$test

path = "/work22/home/redahiro/analysis/OASIS/sceQTL/reQTL/"

##------------------------------##
##        Import data           ##
##------------------------------##


###--- load vcf files ---###

# chrX: exclude autosomal
geno_chrX <- fread(paste0(path,"PseudoBulk/tensorQTL/",cluster,"/",module,"/eSNP/vcf/chrX_Cate",num_cate,".vcf.gz"), header=TRUE, sep='\t', skip='#CHROM') %>%   # data.frame
          select(-c('#CHROM':POS,REF:FORMAT))
colnames_geno_chrX <- colnames(geno_chrX)
colnames(geno_chrX) <- c("variant_id",colnames_geno_chrX[2:length(colnames_geno_chrX)])

geno_chrX <- geno_chrX %>% column_to_rownames("variant_id") %>% t() %>% 
                 as.data.frame() %>% rownames_to_column("SG_ID") %>% as_tibble 


###--- eGene & eSNP paris ---#
eGene_SNP <- read_delim(paste0(path,"PseudoBulk/tensorQTL/",cluster,"/",module,"/Results/summary_",test,"_",module,"_Cate.",num_cate,".txt"),"\t",col_names=T) %>%
               # filter(FDR<0.05) %>% 
               mutate(chr=str_sub(snp,start=1,end=4)) %>% filter(chr == "chrX") %>%   # only retain chrX
               select(gene, snp)

print(paste0("------------ Number of Gene-SNP test: ", dim(eGene_SNP)[1]," ------------"))

###--- Covariates ---###
# Base phenotype
clinical_info <- read_csv("/work22/home/redahiro/analysis/OASIS/phenotype/Phenotype_COV.info_OASIS_eQTL_cov88_hc146.csv") %>%
                      select(SG_ID:Version)

# PC DNA
PC_g <- read_csv("/work22/home/redahiro/analysis/OASIS/phenotype/PCA_genotype_ASA_ALL_cov88_hc148.csv") %>% select(SG_ID,PC1_g,PC2_g)

# Meta info
meta_list <- c("barcode","SG_ID","nCount_RNA","percent.mt")
meta_info <- read_delim(paste0(path,"Data/scRNAseq_base/",cluster,"_metainfo.txt"), "\t", col_names=T) %>% 
                    mutate(SG_ID = if_else(orig.ident == "SG_CO_OK_00021R", "SG_CO_OK_00021", orig.ident)) %>% select(meta_list)

# Cell state info
cell_state_info <- read_delim(paste0(path,"Data/scRNAseq_base/",cluster,"_cellstate.txt"), "\t", col_names=T)

#--- Cov Merge ---#
cov <- cell_state_info %>% 
             inner_join(meta_info, by="barcode") %>%
             inner_join(PC_g, by ="SG_ID") %>%
             inner_join(clinical_info, by="SG_ID") %>%
             mutate(logUMI=log(nCount_RNA))                     # log-transformed of UMI counts

###--- Expression ---###
exp <- read_delim(paste0(path,"SingleCell/Data/",cluster,"_",module,"_UMIdata_update.tsv.gz"), "\t", col_names=T)


#####--------------------------------------------------#####
#####         NB single-cell eQTL analysis             #####
#####--------------------------------------------------#####

lrt_df <- list()
tbl_out_null <- list()
tbl_out_full <- list()

for(i in c(1:dim(eGene_SNP)[1])){
    
    # pickup gene*variant paris 
    GENE = eGene_SNP[i,1] %>% .$gene    
    SNP = eGene_SNP[i,2] %>% .$snp

    print(paste0("Start NB model of ",GENE,": ", Sys.time())) 

    # creating analysis data.frame
    sub_exp <- exp %>% select(barcode,GENE) %>% setnames(c("barcode","E"))
    sub_geno <- geno_chrX %>% select(SG_ID,SNP) %>% setnames(c("SG_ID","G_base")) %>%
                    mutate(G = if_else(G_base %in% c("0","0|0"), "0",
                               if_else(G_base %in% c("1","1|1"), "2", "1"))) %>%
                    select(SG_ID,G)

    data <- sub_exp %>% 
              inner_join(cov, by="barcode") %>%
              inner_join(sub_geno, by="SG_ID") %>%
              mutate_at(vars(G),as.numeric) %>%
              as.data.frame() 
    
    ###--- cell-state*G interaction analysis ---###
    # setting model
    null_form <- paste0("E ~ G + Status + Age + Sex + Version + logUMI + percent.mt + PC1_g + PC2_g + (1 | SG_ID) + ", 
                         paste0("PC_",seq(1,15),collapse=" + "), "+", paste0("harmony_",seq(1,15),collapse=" + "))

    full_form <- paste0("E ~ G + Status + Age + Sex + Version + logUMI + percent.mt + PC1_g + PC2_g + (1 | SG_ID) + ",
                         paste0("PC_",seq(1,15),collapse=" + "), "+", paste0("harmony_",seq(1,15),collapse=" + "), "+", paste0("G*harmony_",seq(1,15),collapse=" + "))
    

    # NB model (about 5min per gene at 180,000 cells)
    null_model <- lme4::glmer.nb(formula = null_form, 
                          nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap"))
   
    full_model <- lme4::glmer.nb(formula = full_form, 
                          nAGQ = 0, data= data, control = glmerControl(optimizer = "nloptwrap"))

    # likelihood ratio test
    model_lrt <- as.data.frame(anova(null_model, full_model))
    lrt_df[[GENE]] <- tibble(Gene=GENE, Snp=SNP, MODEL=row.names(model_lrt), as_tibble(model_lrt))
    
    # Wald tests and effect size information
    out_null <- summary(null_model)$coefficients
    colnames(out_null)<-c("BETA","SE","Z","P")
    tbl_out_null[[GENE]] <- tibble(Gene=GENE, Snp=SNP, Variable=row.names(out_null),as_tibble(out_null))

    out_full <- summary(full_model)$coefficients
    colnames(out_full)<-c("BETA","SE","Z","P")
    tbl_out_full[[GENE]] <- tibble(Gene=GENE, Snp=SNP, Variable=row.names(out_full),as_tibble(out_full))
  
    print(paste0("Finished NB model of ",GENE,": ", Sys.time()))
  }

print(paste0("------------------- ALL DONE: ", Sys.time(), " -------------------"))

LRT <- do.call(rbind, lrt_df)
OUT_NULL <- do.call(rbind, tbl_out_null)
OUT_FULL <- do.call(rbind, tbl_out_full)

# Export
path_out = "/work22/home/redahiro/analysis/OASIS/sceQTL/reQTL/SingleCell/Results/"

write_tsv(LRT,paste0(path_out,cluster,"_",module,"_",test,"_ChrX_int_assoc_NB.LRT.tsv.gz"))
write_tsv(OUT_FULL,paste0(path_out,cluster,"_",module,"_",test,"_ChrX_int_assoc_NB_full.summary.tsv.gz"))
write_tsv(OUT_NULL,paste0(path_out,cluster,"_",module,"_",test,"_ChrX_int_assoc_NB_null.summary.tsv.gz"))
