#--------------------------------------------------------------------------------------------#
# Colocalization using coloc: cis-eQTL
# env: GWAS
# dir: /work22/home/redahiro/analysis/OASIS/GWAS
# date: 2023.11.09
#--------------------------------------------------------------------------------------------#

# load library
library(tidyverse)
library(coloc)
library(optparse)

# options parser
option_list = list(
  make_option("--cell_type", type="character", default=FALSE,
              help="Major Six cell types: CD4T or CD8T or NK or B or MD."),
  make_option("--cell_level", type="character", default=FALSE,
              help="Manual annotation: l1 or l2 or l3."),
  make_option("--pheno", type="character", default=FALSE,
              help="Pheno type of GWAS."))
opt = parse_args(OptionParser(option_list=option_list))

# base_setting
cell_type <- opt$cell_type                 
cell_level <- opt$cell_level
PHENO <- opt$pheno

path="/work22/home/redahiro/analysis/OASIS/GWAS/"
path2="/work22/home/redahiro/analysis/OASIS/sceQTL/PseudoBulk/"
coloc_abf_names <- c("snps","H0","H1","H2","H3","H4")

# lead_variant
lead_variant <- read_csv(paste0("coloc/lead_variant/",PHENO,"_leadvariant_for_coloc.csv")) %>%
                      select(CHR,POS,rsID)

# OASIS cluster
target_cell_list <- read_csv("/work22/home/redahiro/analysis/OASIS/sceQTL/l1_l2_l3_cluster.csv") %>%
                       filter(l1 == cell_type) %>%
                       select(cell_level) %>% setNames("Cluster") %>% .$Cluster %>% unique

sub_summary <- list()
num_variant <- lead_variant$rsID %>% unique %>% length

#-----------------------------------------------------------#
#                      Import base data                     #
#-----------------------------------------------------------#
GWAS_base <- read_delim(paste0("/work22/home/redahiro/analysis/OASIS/GWAS/",PHENO,"/sumstats_for_coloc.txt"),"\t",col_names=T)

#----------------------------------------------------------#
#                       per cell type                      #
#----------------------------------------------------------#
for(t_cell in target_cell_list){
print(paste0("#####--- TARGET CLUSTER:", t_cell, "---#####"))
SampleN_0 <- read_delim(paste0(path2,cell_type,"/",cell_level,"/Covariates/",t_cell,"_Gene.1_Cell.10_no.flip_Cov.txt"),"\t",col_names=T) %>% colnames() %>% length
SampleN <- SampleN_0-1

eQTL_base <- read_delim(paste0(path2,cell_type,"/",cell_level,"/nominal/",t_cell,"_PC15_MAF0.05.cis_nominal.txt.gz"),"\t",col_names=T) %>%
            separate(variant_id,into=c("CHR_0","POS","REF","ALT"),sep="_") %>%
            mutate(CHR=str_sub(CHR_0,start=4,end=-1),
            	     CHR=if_else(CHR=="X","23",CHR),
            	     N=SampleN,                                                    
            	     MAF=af, P=pval_nominal,beta=slope,betase=slope_se) %>%
            filter(!betase %in% NA) %>%
            mutate_at(vars(CHR,POS),as.numeric) %>%
            select(CHR,POS,REF,ALT,beta,betase,P,N,MAF,gene)                    
#-----------------------------------------------------------#
#                 for colocalization anlaysis               #
#-----------------------------------------------------------#
for(i in c(1:num_variant)){
print(paste0("#####--- START:", i," variant"))

# setting target variant
set <- lead_variant[i,]
pheno=PHENO
rsID=set$rsID
chr=set$CHR
pos=set$POS

# filtering Target regions
eQTL_r <- eQTL_base %>% filter(CHR==chr) %>% filter(POS>pos-500000 & POS<pos+500000) %>% mutate(snpID = paste0(CHR,"_",POS,"_",REF,"_",ALT)) %>% mutate(varbeta = betase*betase)
GWAS <- GWAS_base %>% filter(CHR==chr) %>% filter(POS>pos-500000 & POS<pos+500000) %>% mutate(snpID = paste0(CHR,"_",POS,"_",ref,"_",alt)) %>% mutate(varbeta = betase*betase)

gene_list <- eQTL_r$gene %>% unique
print(paste0("   number of target genes: ", length(gene_list)))   
  for(j in c(1:length(gene_list))){
    GENE <- gene_list[j]
    eQTL <- eQTL_r %>% filter(gene==GENE) %>% select(-gene)
    
    # check: CHR & POS & REF & ALT
    overlap <- GWAS %>% inner_join(eQTL,by="POS",suffix=c("_g","_e")) %>%
                  filter(ref==REF & alt==ALT | ref==ALT & alt==REF) %>%
                  mutate(beta = if_else(ref==REF, beta_e, -beta_e)) %>%
                  distinct(snpID_g, .keep_all=TRUE) %>%
                  distinct(snpID_e, .keep_all=TRUE) 
   
    list_e <- list()
    list_e[["beta"]] <- overlap$beta
    list_e[["varbeta"]] <- overlap$varbeta_e
    list_e[["snp"]] <- overlap$snpID_g
    list_e[["position"]] <- overlap$POS
    list_e[["type"]] <- "quant"
    list_e[["N"]] <- SampleN
    list_e[["MAF"]] <- overlap$MAF
    
    list_g <- list()
    list_g[["beta"]] <- overlap$beta_g
    list_g[["varbeta"]] <- overlap$varbeta_g
    list_g[["snp"]] <- overlap$snpID_g
    list_g[["position"]] <- overlap$POS
    list_g[["type"]] <- "cc"
    
    overlap_snp <- dim(overlap)[1]
    if(overlap_snp>0){
           res <- coloc.abf(dataset1=list_e, dataset2=list_g)    
           sub_summary[[paste0(t_cell,"_",rsID,"_",GENE)]] <- res$summary %>% as.data.frame %>% t() %>% as.data.frame %>% as_tibble %>% setNames(coloc_abf_names) %>% 
                                     mutate(Cell=t_cell,Pheno=pheno,rsID=rsID,CHR=chr,POS=pos,gene=GENE) %>% select(Cell,Pheno:gene,everything())
      }
    }
  }
  rm(eQTL_base)
}

#-----------------------------------------------------------#
#                        Export                             #
#-----------------------------------------------------------#
summary_all <- do.call(rbind,sub_summary)
summary_all <- summary_all %>% filter(snps>100)
write.table(summary_all, paste0(path,"coloc/results/cis_nominal/",cell_level,"/", PHENO,"_",cell_type,"_summary_all.txt"), quote=F, sep="\t", row.names=F, col.names=T)

# pick-up top variant of each variants
top_sub <- list()
for(i in target_cell_list){
for(j in c(lead_variant$rsID)){
  sub <- summary_all %>% filter(Cell==i & rsID==j) %>% arrange(-H4)
  top_sub[[paste0(i,"_",j)]] <- sub[1,]
 }
}

summary_top <- do.call(rbind,top_sub) %>% arrange(-H4)
write.table(summary_top, paste0(path,"coloc/results/cis_nominal/",cell_level,"/", PHENO,"_",cell_type,"_summary_top.txt"), quote=F, sep="\t", row.names=F, col.names=T)
