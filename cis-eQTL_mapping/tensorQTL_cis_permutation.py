#------------------------------------------------------------------------------------------------#
# tensorQTL cis permutation        
# env: tensorQTL
# dir: /work22/home/redahiro/analysis/OASIS/sceQTL/PseudoBulk/
# date: 2023.09.13
#------------------------------------------------------------------------------------------------#

import pandas as pd
import numpy as np
import tensorqtl
from tensorqtl import genotypeio, cis, trans
from optparse import OptionParser
import time as tm

#--------------------------------------------------#
#                   optparse                       #
#--------------------------------------------------#
parser = OptionParser()
parser.add_option("-l","--l1_cluster", dest="l1_cluster")                      # major cell type
parser.add_option("-r","--cluster_resolution", dest="cluster_resolution")      # L1, L2, LOneK1K
parser.add_option("-p","--pc_num",
    default="15",type="int",dest="pc_num")
parser.add_option("-m","--maf",
    default="0.05",type="float",dest="maf")

(options,args) = parser.parse_args()

l1_cluster = options.l1_cluster
cluster_resolution = options.cluster_resolution
pc_num = options.pc_num
maf = options.maf

print(l1_cluster)
print(cluster_resolution)
print(pc_num)
print(maf)

## target cluster setting
table_base = pd.read_csv("/work22/home/redahiro/analysis/OASIS/sceQTL/l1_l2_l3_cluster.csv", sep=',')
table = table_base.loc[table_base.l1==l1_cluster,:]
table = table.loc[:,cluster_resolution] 
target_cluster = table.unique()

print(f'target cluster: {target_cluster}')

## Making path
base_path = "/work22/home/redahiro/analysis/OASIS/sceQTL/PseudoBulk/"
path = f'{base_path}{l1_cluster}/{cluster_resolution}/'

# chr list: only autosomal
chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
	        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
            'chr21','chr22']
KF_ID = {'SG_XXXX'}           # sample ID removed from chrX analysis

#--------------------------------------------------#
#          eQTL mapping: permutation               #
#--------------------------------------------------#

for i in target_cluster:

    print("starting {0}, {1}".format(i, tm.ctime()))
    
	##--- Phenotype & Covariates ---##
    phenotype_bed_file = f'{path}BED/{i}_Gene.1_Cell.10_Input_no.flip_ALL.bed.gz'
    covariates_file = f'{path}Covariates/{i}_Gene.1_Cell.10_no.flip_Cov.txt'             
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T                # samples x covariates
    if len(covariates_df.index.intersection(KF_ID)) == 1 :                               # dimension fixed for phenotype of chrX
        covariates_df_X = covariates_df.drop('SG_XXXX', axis=0)
    else :
        covariates_df_X = covariates_df                  


	##--- Genotypeds ---##

    # Autosomal
    plink_prefix_path = '/work22/home/redahiro/analysis/OASIS/WGS/WGS_202308/plink/Merged_Auto_WGS_202308'
    pr = genotypeio.PlinkReader(plink_prefix_path)
    # Chr X
    plink_prefix_path_X = '/work22/home/redahiro/analysis/OASIS/WGS/WGS_202308/plink/X_nonPAR.S4phased'
    pr_X = genotypeio.PlinkReader(plink_prefix_path_X)

	## load genotypes and variants into data frames
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    genotype_df_X = pr_X.load_genotypes()
    variant_df_X = pr_X.bim.set_index('snp')[['chrom', 'pos']]

	##--- covariate setting ---##
    cov_list = ('Status','Age','Sex','Version','PC1_g','PC2_g',
                'PC1_r','PC2_r','PC3_r','PC4_r','PC5_r','PC6_r','PC7_r','PC8_r','PC9_r','PC10_r',
                'PC11_r','PC12_r','PC13_r','PC14_r','PC15_r','PC16_r','PC17_r','PC18_r','PC19_r','PC20_r',
                'PC21_r','PC22_r','PC23_r','PC24_r','PC25_r','PC26_r','PC27_r','PC28_r','PC29_r','PC30_r')
    
    # selecting PC_rna number
    cov_list = cov_list[0:pc_num+6]
    covariates_df = covariates_df.loc[:, cov_list]
    covariates_df_X = covariates_df_X.loc[:, cov_list]

    print(f'list of covariates: {covariates_df.columns}')

    ##--- phenotype df needs to be subsetted ---##
    phenotype_df_sub = phenotype_df.loc[phenotype_pos_df['chr'].isin(chr_list),:]
    phenotype_pos_df_sub = phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(chr_list),:]

    phenotype_df_X = phenotype_df.loc[phenotype_pos_df['chr'].isin(['chrX']),:]
    if len(covariates_df.index.intersection(KF_ID)) == 1 :                        
        phenotype_df_X = phenotype_df_X.drop('SG_CO_OK_00057', axis=1)
    else :
        phenotype_df_X = phenotype_df_X        
    phenotype_pos_df_X = phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(["chrX"]),:]

    ##-----------------------------------##
    ##            permutation            ##
    ##-----------------------------------##
    prefix = f'permu/{i}_PC{pc_num}_MAF{maf}'
    # Autosomal 
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df_sub, phenotype_pos_df_sub, covariates_df, maf_threshold=maf)
    # chrX
    cis_df_X = cis.map_cis(genotype_df_X, variant_df_X, phenotype_df_X, phenotype_pos_df_X, covariates_df_X, maf_threshold=maf)

    # Merge
    cis_all = pd.concat([cis_df,cis_df_X],axis=0)

    # calculate genome-wide qvalues  
    tensorqtl.calculate_qvalues(cis_all)

    # eGene check
    eGene_n = len(cis_all.query('qval<0.05'))
    print(f'eGene number: {eGene_n}')
    list_sub = [i,eGene_n]
    eGene_n = np.array(list_sub)
    eGene_df = pd.DataFrame(eGene_n)

    # gene name added
    gene_id = pd.read_csv("/work22/home/redahiro/reference_data/scRNAseq/making_bed/features.tsv.gz",
	                      sep='\t', header=None, names=('phenotype_id','gene','GE'))
    gene_id = gene_id.loc[:,['phenotype_id','gene']]
    cis_all = cis_all.merge(gene_id, on=['phenotype_id'], how='left')

    # reindex of columns
    cis_all = cis_all.reindex(columns=['phenotype_id','gene','variant_id','tss_distance','qval','pval_nominal','slope','slope_se','af',
	                                 'pval_beta','pval_perm','pval_nominal_threshold','num_var','beta_shape1','beta_shape2','true_df','pval_true_df','ma_samples','ma_count'])

	# export
    cis_all.to_csv(f'{path}{prefix}_Cell.10_top_assoc_chr1_23.txt.gz', sep='\t', index=False, compression = "gzip")
    eGene_df.to_csv(f'{path}{prefix}_Cell.10_eGene_number_chr1_23.txt', sep='\t', index=False)

    print("done {0}, {1}".format(i, tm.ctime()))

print('cis-QTL Done')
