#------------------------------------------------------------------------------------------------#
# Roadmap: Functional enrichment analysis of cis-eSNPs    
# env: Scanpy
# dir: /work22/home/redahiro/analysis/OASIS/sceQTL/PseudoBulk/anno/
# date: 2023.09.11
#------------------------------------------------------------------------------------------------#

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

l1_cluster = options.l1_cluster
cluster_resolution = options.cluster_resolution

print(l1_cluster)
print(cluster_resolution)

## target cluster setting
table_base = pd.read_csv("/work22/home/redahiro/analysis/OASIS/sceQTL/l1_l2_l3_cluster.csv", sep=',')
table = table_base.loc[table_base.l1==l1_cluster,:]
table = table.loc[:,cluster_resolution] 
target_cluster = table.unique()

print(f'target cluster: {target_cluster}')

## Making path
base_path = '/work22/home/redahiro/analysis/OASIS/sceQTL/PseudoBulk/'
path_load = f'{base_path}{l1_cluster}/{cluster_resolution}/'
path_out = f'{base_path}anno/{cluster_resolution}/'
path_anno = '/work22/home/redahiro/reference_data/epigenomics/Roadmap/hg38_anno'

## pickup roadmap celltype list
roadmap_ID_PBMC = ['E029','E032','E046','E038','E037','E044','E047','E048']
print('Prepareation DONE')

#---------------------------------------------------#
#    Enrichment anlaysis of Roadmap data (hg38)     #
#---------------------------------------------------#

print('Enrichment analysis of PBMC: START')

# null dataframe for summary data
df0 = []

for cluster in target_cluster:
    print(f'{cluster} START')
    
    # Import permutation & nominal eQTL results
    eSNP = pd.read_csv(f'{path_load}permu/{cluster}_PC15_MAF0.05_Cell.10_top_assoc_chr1_23.txt.gz', sep="\t")
    eSNP = eSNP.query("qval<0.05")                      # only significant-eSNPs
    bg = pd.read_csv(f'{path_load}nominal/{cluster}_PC15_MAF0.05.cis_nominal.txt.gz', sep="\t")
    
    # Making bed files for analysis
    eSNP["chr"] = eSNP.variant_id.str.split("_").str[0]
    eSNP["start"] = eSNP.variant_id.str.split("_").str[1].astype(int) -1
    eSNP["end"] = eSNP.start + 1
    variants_bed = eSNP.iloc[:,-3:]
    variants_bed.to_csv(f'{path_out}/tmp_eSNP/{cluster}_eSNP_positions.bed.gz', sep='\t', index=False, header=False, compression = "gzip")

    bg["chr"] = bg.variant_id.str.split("_").str[0]
    bg["start"] = bg.variant_id.str.split("_").str[1].astype(int) -1
    bg["end"] = bg.start + 1
    variants_bed_bg = bg.iloc[:,-3:]
    variants_bed_bg.to_csv(f'{path_out}/tmp_eSNP/{cluster}_bg_positions.bed.gz', sep='\t', index=False, header=False, compression = "gzip")

    fn_e = f'{path_out}tmp_eSNP/{cluster}_eSNP_positions.bed.gz'
    ebed = pybedtools.BedTool(fn_e)

    fn_bg = f'{path_out}tmp_eSNP/{cluster}_bg_positions.bed.gz'
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
        case_n = eSNP.shape[0]
        case_p = len(it_pro)
        case_non_p = case_n - case_p
        case_e = len(it_enh)
        case_non_e = case_n - case_e
        
        bg_n = bg.shape[0]
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
    	       'L1_Cluster': [l1_cluster],
               'Resolution': [cluster_resolution],
               'Cluster': [cluster], 
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
    print(f'{cluster} Done')
df0 = pd.concat(df0)
df0.to_csv(f'{path_out}/tmp_results_eSNP/{l1_cluster}_results.txt', sep='\t', index=False, header=True)
print('Enrichment analysis of PBMC: DONE')


###--- ALL Roadmap data ---###
fns = glob.glob('/work22/home/redahiro/reference_data/epigenomics/Roadmap/hg38/*.bed.gz')
roadmap_ID_ALL = pd.Series(fns).apply(lambda x: x.replace("/work22/home/redahiro/reference_data/epigenomics/Roadmap/hg38/", "").replace("_18_core_K27ac_hg38lift_mnemonics.bed.gz",""))

print('Enrichment analysis of ALL Data: START')

# null dataframe for summary data
df1 = []

for cluster in target_cluster:
    print(f'{cluster} START')
    
    # Import permutation & nominal eQTL results
    fn_e = f'{path_out}tmp_eSNP/{cluster}_eSNP_positions.bed.gz'
    ebed = pybedtools.BedTool(fn_e)

    fn_bg = f'{path_out}tmp_eSNP/{cluster}_bg_positions.bed.gz'
    bgbed = pybedtools.BedTool(fn_bg)

    #--- Enrichment analysis for each celltype ---# 
    for i in roadmap_ID_ALL:
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
        case_n = eSNP.shape[0]
        case_p = len(it_pro)
        case_non_p = case_n - case_p
        case_e = len(it_enh)
        case_non_e = case_n - case_e

        bg_n = bg.shape[0]
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
               'L1_Cluster': [l1_cluster],
               'Resolution': [cluster_resolution],
               'Cluster': [cluster], 
               'RoadMap': [i],
               'Odds_p':[fisher_p[0]],
               'SE_sqrt_p':[se_sqrt_p],
               'Pvalue_p':[fisher_p[1]],
               'Odds_e':[fisher_e[0]],
               'SE_sqrt_e':[se_sqrt_e],
               'Pvalue_e':[fisher_e[1]],
            })

        df1.append(df)
        print(f'    {i} Done')
    print(f'{cluster} Done')
df1 = pd.concat(df1)
df1.to_csv(f'{path_out}/tmp_results_eSNP_ALL/{l1_cluster}_results.txt', sep='\t', index=False, header=True)
print('Enrichment analysis of ALL Data: DONE')
