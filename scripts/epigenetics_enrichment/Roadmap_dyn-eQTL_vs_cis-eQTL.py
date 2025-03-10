#----------------------------------------------------------------------------------------------#
# Roadmap: Functional enrichment analysis of eSNP
# env: Scanpy
# dir: /work22/home/redahiro/analysis/OASIS/sceQTL/reQTL/
# date: 2023.09.27 → update 2024.02.16 (including Sex as covariates & adding L3 cluster)
#----------------------------------------------------------------------------------------------#

import pandas as pd
import numpy as np
import time as tm
from scipy import stats
from optparse import OptionParser
import glob
import pybedtools
import scipy.stats as st

# optparse
parser = OptionParser()
parser.add_option("-l","--l1_cluster", dest="l1_cluster")
parser.add_option("-r","--cluster_resolution", dest="cluster_resolution")

(options,args) = parser.parse_args()

base_cluster = options.l1_cluster
base_cluster = "Mono"
target_cluster = "MD"
module = "HM_IR"   # "HM_IR" "GO_IFNG" "TYPE1_ISG"
num_cate = 10

## Making path
base_path = '/work22/home/redahiro/analysis/OASIS/sceQTL/reQTL/PseudoBulk/tensorQTL/'
path_base_load = f'/work22/home/redahiro/analysis/OASIS/sceQTL/PseudoBulk/{base_cluster}/l1/permu/'
path_target_load = f'{base_path}{target_cluster}/{module}/Results/'
path_out = f'{base_path}anno/BED/'
path_out2= f'{base_path}anno/Summary'
path_anno = '/work22/home/redahiro/reference_data/epigenomics/Roadmap/hg38_anno'

## pickup roadmap celltype list
roadmap_ID_PBMC = ['E029','E032','E046','E038','E037','E044','E047','E048']
test_list = ['linear','quadratic']

print('Prepareation DONE')

#---------------------------------------------------#
#    Enrichment anlaysis of Roadmap data            #
#---------------------------------------------------#

print('Enrichment analysis of PBMC: START')

for test in test_list:
    # null dataframe for summary data
    df0 = []
    # Import permutation & nominal eQTL results
    deSNP = pd.read_csv(f'{path_target_load}summary_{test}_{module}_Cate.{num_cate}_update.txt', sep="\t") 
    deSNP = deSNP.query("FDR<0.05")          # pick-up significant dynamic-eQTLs
    eSNP_bg = pd.read_csv(f'{path_base_load}{base_cluster}_PC15_MAF0.05_Cell.10_top_assoc_chr1_23.txt.gz', sep="\t")
    eSNP_bg = eSNP_bg.query("qval<0.05")     # pick-up significant cis-eQTLs
   
    # Making bed files for analysis
    deSNP["chr"] = deSNP.snp.str.split("_").str[0]
    deSNP["start"] = deSNP.snp.str.split("_").str[1].astype(int) -1
    deSNP["end"] = deSNP.start + 1
    variants_bed = deSNP.iloc[:,-3:]
    variants_bed.to_csv(f'{path_out}/dyneSNP_{module}_{test}_{num_cate}_positions.bed', sep='\t', index=False, header=False)

    eSNP_bg["chr"] = eSNP_bg.variant_id.str.split("_").str[0]
    eSNP_bg["start"] = eSNP_bg.variant_id.str.split("_").str[1].astype(int) -1
    eSNP_bg["end"] = eSNP_bg.start + 1
    variants_bed_bg = eSNP_bg.iloc[:,-3:]
    variants_bed_bg.to_csv(f'{path_out}/bg_positions.bed', sep='\t', index=False, header=False)

    fn_e = f'{path_out}/dyneSNP_{module}_{test}_{num_cate}_positions.bed'
    ebed = pybedtools.BedTool(fn_e)  
    fn_bg = f'{path_out}/bg_positions.bed'
    bgbed = pybedtools.BedTool(fn_bg)
    
    #--- Enrichment analysis for each celltype ---# 
    for i in roadmap_ID_PBMC:
        print(f'    {i} start')
        # promoter
        RD_pro = pybedtools.example_bedtool(f'{path_anno}/{i}_promoter.bed')
        it_pro = ebed.intersect(RD_pro)
        itbg_pro = bgbed.intersect(RD_pro)
        # enhancer
        RD_enh = pybedtools.example_bedtool(f'{path_anno}/{i}_enhancer.bed')
        it_enh = ebed.intersect(RD_enh)
        itbg_enh = bgbed.intersect(RD_enh)
        # calculate OR
        case_n = deSNP.shape[0]
        case_p = len(it_pro)
        case_non_p = case_n - case_p
        case_e = len(it_enh)
        case_non_e = case_n - case_e
        
        bg_n = eSNP_bg.shape[0]
        bg_p = len(itbg_pro)
        bg_non_p = bg_n - bg_p
        bg_e = len(itbg_enh)
        bg_non_e = bg_n - bg_e

        table_p = [[case_p,case_non_p],[bg_p,bg_non_p]]
        fisher_p = st.fisher_exact(table_p, alternative='two-sided')
        table_e = [[case_e,case_non_e],[bg_e,bg_non_e]]
        fisher_e = st.fisher_exact(table_e, alternative='two-sided')
        se_sqrt_p = np.sqrt(1/case_p + 1/case_non_p + 1/bg_p + 1/bg_non_p)
        se_sqrt_e = np.sqrt(1/case_e + 1/case_non_e + 1/bg_e + 1/bg_non_e)

        # Making df
        df = pd.DataFrame({
               'Cluster': [target_cluster],
    	       'Module': [module],
    	       'RoadMap': [i],
    	       'Odds_p':[fisher_p[0]],
    	       'SE_sqrt_p':[se_sqrt_p],
               'Pvalue_p':[fisher_p[1]],
    	       'Odds_e':[fisher_e[0]],
    	       'SE_sqrt_e':[se_sqrt_e],
               'Pvalue_e':[fisher_e[1]],
    	    })        
        df0.append(df)
        print(f'    {i} Done')
    df0 = pd.concat(df0)
    df0.to_csv(f'{path_out2}/{target_cluster}_{module}_{test}_{num_cate}_results_update.txt', sep='\t', index=False, header=True)   
    print(f'{test} Done')

print('Enrichment analysis of PBMC: DONE')

